/*
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

#include <cassert>
#include <cmath>

#include "PowerLawNu.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam::heatTransferModels
{
defineTypeNameAndDebug(PowerLawNu, 0);
addToRunTimeSelectionTable(heatTransferModel, PowerLawNu, dictionary);
} // namespace Foam::heatTransferModels

Foam::heatTransferModels::PowerLawNu::PowerLawNu(
    const dictionary &dict, const phaseInterface &interface,
    const bool registerObject)
    : heatTransferModel{dict, interface, registerObject},
      interface_(
          interface.modelCast<heatTransferModel, dispersedPhaseInterface>()),
      coeff_a_{"a", dimless, dict}, coeff_b_{"b", dimless, dict},
      coeff_c_{"c", dimless, dict},
      superficial_Re_{dict.lookupOrDefault("superficialRe", false)},
      residual_Re_{dict.lookupOrDefault("residualRe", 1e-6)}
{
    if (coeff_a_.value() <= 0) {
        WarningInFunction << "Factor a = " << coeff_a_.value()
                          << " seems non-physical for PowerLawNu" << endl;
    }
    if (coeff_b_.value() < 0) {
        WarningInFunction << "Factor b = " << coeff_b_.value()
                          << " seems non-physical for PowerLawNu" << endl;
    }
    if (coeff_c_.value() < 0) {
        WarningInFunction << "Factor c = " << coeff_c_.value()
                          << " seems non-physical for PowerLawNu" << endl;
    }
    if (residual_Re_ <= 0) {
        FatalErrorInFunction << "Residual Re must be positive, found "
                             << residual_Re_ << exit(FatalError);
    }
}

Foam::tmp<Foam::volScalarField> Foam::heatTransferModels::PowerLawNu::K(
    const scalar residualAlpha) const
{
    auto Re = max(interface_.Re(), residual_Re_);

    volScalarField alpha_bounded(max(interface_.dispersed(), residualAlpha));
    if (superficial_Re_) {
        Re.ref() *= alpha_bounded;
    }

    volScalarField Nu(coeff_a_ * pow(Re, coeff_b_) *
                      pow(interface_.Pr(), coeff_c_));
    assert(min(Nu).value() > 0);

    return 6 * alpha_bounded * interface_.continuous().thermo().kappa() * Nu /
        sqr(interface_.dispersed().d());
}