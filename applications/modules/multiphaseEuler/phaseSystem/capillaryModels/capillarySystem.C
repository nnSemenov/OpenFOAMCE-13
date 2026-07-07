/*---------------------------------------------------------------------------*\
License
    This file is part of Mikeno. Mikeno (https://github.com/nnSemenov/Mikeno) is
    a fork of OpenFOAM (https://github.com/OpenFOAM/OpenFOAM-dev).

    Mikeno is free software: you can redistribute it and/or modify it  under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    Mikeno is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License
    along with Mikeno.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include <Time.H>
#include <fluidThermo.H>
// #include <fvcGrad.H>
// #include <fvcSnGrad.H>
#include <fvc.H>

#include "capillarySystem.H"
#include "correctFixedFluxBCs.H"

#include <cassert>
#include <filesystem>

namespace Foam
{
defineTypeNameAndDebug(capillarySystem, 0)
} // namespace Foam

using namespace Foam;

capillarySystem::capillarySystem(const phaseSystem &phases)
    : IOdictionary(IOobject{dictName, Foam::Time::constant(), phases.local(),
                            phases.mesh(), IOobject::MUST_READ,
                            IOobject::NO_WRITE, true}),
      phase_system_{phases},
      reference_phase_{this->lookup<word>("referencePhase")}
{
    if (not phases.found(reference_phase_)) {
        FatalIOErrorInFunction(*this)
            << "Reference phase " << reference_phase_
            << " is missing from phase system" << exit(FatalIOError);
    }

    read_capillary_models(*this);
}

std::vector<const phaseModel *> capillarySystem::fluid_phases() const
{
    auto &phases = phase_system_.phases();
    std::vector<const phaseModel *> ret;
    ret.reserve(phases.size());
    for (auto &phase : phases) {
        if (phase.stationary()) {
            continue;
        }

        auto &thermo = phase.thermo();
        if (dynamic_cast<const fluidThermo *>(&thermo) == nullptr) {
            // Thermo is not fluid thermo, considered to be solid
            continue;
        }
        ret.emplace_back(&phase);
    }
    return ret;
}

std::vector<word> capillarySystem::fluid_phase_names() const
{
    auto phases = fluid_phases();
    std::vector<word> ret;
    ret.reserve(phases.size());
    for (auto model : phases) {
        ret.emplace_back(model->name());
    }
    return ret;
}

void capillarySystem::read_capillary_models(const dictionary &dict) &
{
    auto phases_todo = fluid_phases();

    const phaseModel &reference = phase_system_.phases()[reference_phase_];

    const dictionary &model_dicts = dict.subDict("capillaryPressure");
    for (const phaseModel *phase : phases_todo) {
        const word phase_name = phase->name();
        if (phase_name == reference_phase_) {
            continue;
        }

        const word key = phase_pair_name(reference_phase_, phase_name);
        const phaseInterface interface{reference, *phase};

        auto &capillary_dict = model_dicts.subDict(key);

        auto model = capillaryPressureModels::capillaryPressureModel::New(
            capillary_dict, interface);
        this->capillaryPressureModels_.emplace(phase_name, std::move(model));
    }
}

volScalarField capillarySystem::pressure_difference_from_common_p(
    const word &phase_i) const
{

    auto &mesh = phase_system_.mesh();
    volScalarField diff_p(IOobject("p-p_" + phase_i, mesh, IOobject::NO_READ,
                                   IOobject::NO_WRITE, false),
                          mesh, dimensionedScalar{dimPressure, scalar{0}});

    if (phase_i == reference_phase_) {
        return diff_p;
    }
    diff_p = this->capillaryPressureModels_.at(phase_i)->capillary_pressure();
    return diff_p;
}

void capillarySystem::warn_if_movable_solids(const word &phase) const
{
    auto &pm = phase_system_.phases()[phase];

    if (phase == reference_phase_) {
        return;
    }

    if (capillaryPressureModels_.find(phase) !=
        capillaryPressureModels_.end()) {
        // considered to be fluid
        return;
    }

    if (pm.stationary()) {
        return;
    }

    WarningInFunction << "Capillary system found movable solid phase " << phase
                      << ", currently not able to process capillary force on "
                         "movable solids correctly."
                      << endl;
}

tmp<volVectorField> capillarySystem::F(const word &phase) const
{
    const word name = "F_cap_" + phase;
    auto ret = fvc::reconstruct(Ff(phase));
    ret.ref().rename(name);
    return ret;
}

tmp<surfaceScalarField> capillarySystem::Ff(const word &phase) const
{

    const volScalarField &alpha = phase_system_.phases()[phase];
    const word name = "Ff_cap_" + phase;
    auto &mesh = phase_system_.mesh();

    warn_if_movable_solids(phase);
    if (phase == reference_phase_) {
        return surfaceScalarField ::New(
            name, mesh,
            dimensionedScalar{dimForce / dimVolume * dimArea, Foam::Zero});
    }

    auto it = capillaryPressureModels_.find(phase);
    if (it == capillaryPressureModels_.end()) {
        // Solid phase. Currently only supports stationary solids
        return surfaceScalarField ::New(
            name, mesh,
            dimensionedScalar{dimForce / dimVolume * dimArea, Foam::Zero});
    }

    auto &cap_model = it->second;
    auto p_sub_pi = cap_model->capillary_pressure();
    auto force = fvc::interpolate(alpha) * fvc::snGrad(p_sub_pi, "snGrad(Pc)") *
        mesh.magSf();

    auto ret = correctFixedFluxBCs(cap_model->interface(), force);
    ret.ref().rename(name);

    return ret;
}

void capillarySystem::add_to_Fs(PtrList<surfaceScalarField> &Fs) const
{
    add_to_Ffs(Fs);
}

void capillarySystem::add_to_Ffs(PtrList<surfaceScalarField> &Ffs) const
{

    assert(Ffs.size() == phase_system_.phases().size());
    //    PtrList<surfaceScalarField> ret{phase_system_.size()};

    for (auto &fluid : fluid_phases()) {
        const word &name = fluid->name();
        auto force = Ff(name);

        addField(*fluid, "Ff", force, Ffs);
    }
}

volScalarField capillarySystem::pressure_difference_between_phases(
    const word &phase1, const word &phase2) const
{
    auto &mesh = phase_system_.mesh();
    volScalarField diff_p(IOobject("p_" + phase1 + "-p_" + phase2, mesh,
                                   IOobject::NO_READ, IOobject::NO_WRITE,
                                   false),
                          mesh, dimensionedScalar{dimPressure, scalar{0}});
    if (phase1 == phase2) {
        return diff_p;
    }

    if (phase1 == reference_phase_) {
        diff_p +=
            this->capillaryPressureModels_.at(phase2)->capillary_pressure();
        return diff_p;
    }
    if (phase2 == reference_phase_) {
        diff_p -=
            this->capillaryPressureModels_.at(phase1)->capillary_pressure();
        return diff_p;
    }

    diff_p = -this->capillaryPressureModels_.at(phase1)->capillary_pressure() +
        this->capillaryPressureModels_.at(phase2)->capillary_pressure();
    return diff_p;
}


autoPtr<capillarySystem> capillarySystem::try_read(
    const phaseSystem &phaseSystem, const fvMesh &mesh)
{
    auto &runTime = mesh.time();
    const bool exist = [&]() {
        IOobject obj(dictName, runTime.caseConstant(), mesh,
                     IOobject::READ_IF_PRESENT, IOobject::NO_WRITE, false);

        return obj.headerOk();
    }();

    if (exist) {
        auto ret = autoPtr<capillarySystem>{new capillarySystem{phaseSystem}};
        assert(mesh.foundObject<capillarySystem>(dictName));

        return ret;
    }

    Info << "Skip capillary system because " << dictName << " doesn't exist."
         << endl;
    return autoPtr<capillarySystem>{nullptr};
}
