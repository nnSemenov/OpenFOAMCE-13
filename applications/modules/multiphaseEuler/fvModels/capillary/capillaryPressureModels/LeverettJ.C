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


#include "LeverettJ.H"
#include "addToRunTimeSelectionTable.H"
#include "averageDiameter.H"
#include "fvcGrad.H"

#include <cassert>

namespace Foam::capillaryPressureModels
{
defineTypeNameAndDebug(LeverettJ, 0);
addToRunTimeSelectionTable(capillaryPressureModel, LeverettJ, dictionary);
} // namespace Foam::capillaryPressureModels

using namespace Foam;
using namespace Foam::capillaryPressureModels;

LeverettJ::LeverettJ(const dictionary &dict, const phaseInterface &interface)
    : capillaryPressureModel{dict, interface},
      liquidName_{dict.lookup<word>("liquid")},
      solidNames_{dict.lookup<List<word>>("solids")},
      E1_{dict.lookup<scalar>("E1")}, J0_{dict.lookup<scalar>("J0")},
      beta_{dict.lookup<scalar>("beta")}, contactAngle_{std::nullopt}
{
    if (solidNames_.empty()) {
        FatalIOErrorInFunction(dict) << "LeverettJ capillary pressure model "
                                        "requires at least one solid phase."
                                     << exit(FatalIOError);
    }
    if (dict.found("contactAngle")) {
        const auto val = dict.lookup<scalar>("contactAngle");

        if (val < 0 or val >= 180) {
            FatalIOErrorInFunction(dict)
                << "Non-physical contact angle " << val
                << ", should be in range [0,180) (degree)"
                << exit(FatalIOError);
        }
        contactAngle_ = val;
    }

    if (not interface.fluid().phases().found(liquidName_)) {
        FatalIOErrorInFunction(dict)
            << "Liquid phase named " << liquidName_
            << " not exist in phase system" << exit(FatalIOError);
    }

    if (interface.phase1().name() != liquidName_ and
        interface.phase2().name() != liquidName_) {
        FatalIOErrorInFunction(dict) << "Assigned liquid phase " << liquidName_
                                     << " doesn't exist in interface "
                                     << interface.name() << exit(FatalIOError);
    }
}

tmp<volScalarField> LeverettJ::capillary_pressure() const
{
    assert(interface_.fluid().phases().found(liquidName_));
    const auto &liq = interface_.fluid().phases()[liquidName_];
    const auto &alpha_liq = dynamic_cast<const volScalarField &>(liq);

    const auto solid_avg = averageDiameter::sauter(
        solidNames_.size(), [&](label solid_idx) -> const phaseModel & {
            const word solid_name = solidNames_[solid_idx];
            return this->interface_.fluid().phases()[solid_name];
        });
    const volScalarField porosity(
        max(1 - solid_avg.alpha, liq.residualAlpha()));

    const volScalarField sat(
        max(min((alpha_liq / porosity), (1 - liq.residualAlpha().value())),
            liq.residualAlpha()));

    auto J = J0_ + beta_ * log((1 - sat) / sat);

    auto sqrt_porosity_over_K =
        solid_avg.alpha / (porosity * solid_avg.d) * sqrt(E1_);

    tmp<volScalarField> Pc =
        sigma() * sqrt_porosity_over_K * J * cos(contactAngle_.value_or(0));
    Pc->rename("Pc_" + interface_.name());
    return Pc;
}