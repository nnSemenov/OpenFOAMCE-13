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

namespace Foam
{
defineTypeNameAndDebug(capillarySystem, 0);
} // namespace Foam

using namespace Foam;

capillarySystem::capillarySystem(const phaseSystem &phases,
                                 const dictionary &dict)
    : IOdictionary(IOobject{"capillarySystem", Time::constant(), phases.mesh(),
                            IOobject::NO_READ, IOobject::NO_WRITE, true},
                   dict),
      phase_system_{phases}, reference_phase_{dict.lookup("referencePhase")}
{
    if (not phases.found(reference_phase_)) {
        FatalIOErrorInFunction(dict)
            << "Reference phase " << reference_phase_
            << " is missing from phase system" << exit(FatalIOError);
    }

    read_capillary_models(dict);
}

std::vector<const phaseModel *> capillarySystem::fluid_phases() const
{
    auto &phases = phase_system_.phases();
    std::vector<const phaseModel *> ret;
    ret.reserve(phases.size());
    forAll(phases, idx)
    {
        if (phases[idx].stationary()) {
            continue;
        }

        auto &thermo = phases[idx].thermo();
        if (dynamic_cast<const fluidThermo *>(&thermo) == nullptr) {
            // Thermo is not fluid thermo, considered to be fluid
            continue;
        }
        ret.emplace_back(&phases[idx]);
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
    auto phases_todo = fluid_phase_names();
    const auto Nf = static_cast<scalar>(phases_todo.size());

    auto &mesh = phase_system_.mesh();
    volScalarField diff_p(IOobject("p-p_" + phase_i, mesh, IOobject::NO_READ,
                                   IOobject::NO_WRITE, false),
                          mesh, dimensionedScalar{dimPressure, scalar{0}});

    if (phase_i not_eq reference_phase_) {

        diff_p += (Nf - 1) / Nf *
            this->capillaryPressureModels_.at(phase_i)->capillary_pressure();
    }
    for (auto &name_j : phases_todo) {
        if (name_j == phase_i) {
            continue;
        }
        if (name_j == reference_phase_) {
            continue;
        }
        diff_p -=
            this->capillaryPressureModels_.at(name_j)->capillary_pressure() /
            Nf;
    }
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
