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

#include "capillarySystem.H"
#include <Time.H>
#include <fluidThermo.H>
#include <fvc.H>
#include <zeroGradientFvPatchFields.H>
#include <dragModel.H>

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
            // Thermo is not fluid thermo, considered to be fluid
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
    this->capillaryPressureModels_.clear();
    this->use_snGrad_ = dict.lookup<bool>("useSnGrad");

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
    auto phases_todo = fluid_phase_names();

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

volVectorField capillarySystem::capillary_force(const word &name) const
{
    auto &alphaL = phase_system_.phases()[name];

    auto &mesh = phase_system_.mesh();
    const word force_name = "F_cap_" + name;

    volVectorField force_i{IOobject{force_name, mesh, IOobject::NO_READ,
                                    IOobject::NO_WRITE, false},
                           mesh, dimensionedVector{dimForce / dimVolume, Zero}};

    auto it = capillaryPressureModels_.find(name);
    if (it == capillaryPressureModels_.end()) {
        return force_i;
    }

    [[maybe_unused]] auto [Pc, dPc_d_ln_alphaL] =
        it->second->capillary_pressure_with_derivative();

    if (use_snGrad_) {
        surfaceScalarField force_f(fvc::interpolate(dPc_d_ln_alphaL) *
                                   fvc::snGrad(alphaL, "snGrad(alpha)") *
                                   mesh.magSf());
        force_i = fvc::reconstruct(force_f);
    }
    else {
        force_i = dPc_d_ln_alphaL * fvc::grad(alphaL, "grad(alpha)");
    }

    return force_i;
}

volScalarField capillarySystem::alpha_diffusivity(const word &name) const
{
    auto &alpha_ref = phase_system_.phases()[reference_phase_];
    auto &alphaL = phase_system_.phases()[name];
    auto &mesh = phase_system_.mesh();

    volScalarField D_eff{IOobject{"D_eff_" + alphaL.name(), mesh,
                                  IOobject::NO_READ, IOobject::NO_WRITE, false},
                         mesh, dimensionedScalar{dimKinematicViscosity, -small},
                         zeroGradientFvPatchScalarField::typeName};

    auto it = capillaryPressureModels_.find(name);
    if (it == capillaryPressureModels_.end()) {
        return D_eff;
    }

    auto [_, dPc_d_ln_alphaL] =
        it->second->capillary_pressure_with_derivative();
    assert(max(dPc_d_ln_alphaL.ref()).value() <= 0);

    auto &K_model = [&]() -> const blendedDragModel & {
        auto &interface = it->second->interface();
        if (phase_system_.foundInterfacialModel<blendedDragModel>(interface)) {
            return phase_system_.lookupInterfacialModel<blendedDragModel>(
                interface);
        }
        phaseInterface interface_inv{interface.phase2(), interface.phase1()};
        return phase_system_.lookupInterfacialModel<blendedDragModel>(
            interface_inv);
    }();
#warning "TODO: Add K_SL"

    auto K_sum =
        (1 + alphaL / max(alpha_ref, alpha_ref.residualAlpha())) * K_model.K();
    assert(min(K_sum.ref()).value() > 0);

    auto add = sqr(alphaL) * dPc_d_ln_alphaL / K_sum;
    D_eff.internalFieldRef() = add->internalField();

    // This diffusivity should be negative, because capillary force is inverse
    // direction than grad(alpha).
    assert(max(D_eff).value() <= 0);
    return D_eff;
}
volScalarField capillarySystem::most_alpha_diffusivity() const
{
    volScalarField ret = alpha_diffusivity(reference_phase_);
    for (auto &[phase, _] : capillaryPressureModels_) {
        auto D_eff = alpha_diffusivity(phase);
        // Take minimum value because diffusivity should be non-positive
        ret.internalFieldRef() =
            min(ret.internalField(), D_eff.internalField());
        assert(max(ret).value() <= 0);
    }
    return ret;
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

    Info << "Skip capillary system because " << dictName << " doesn't exist.";
    //         << endl;
    return autoPtr<capillarySystem>{nullptr};
}
