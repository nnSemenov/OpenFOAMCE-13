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

#include "Gnielinski.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam::heatTransferModels
{
defineTypeNameAndDebug(Gnielinski, 0);
addToRunTimeSelectionTable(heatTransferModel, Gnielinski, dictionary);
} // namespace Foam::heatTransferModels

Foam::heatTransferModels::Gnielinski::Gnielinski(
    const dictionary &dict, const phaseInterface &interface,
    const bool registerObject)
    : heatTransferModel{dict, interface, registerObject},
      interface_(
          interface.modelCast<heatTransferModel, dispersedPhaseInterface>()),
      sphere_{dict.lookup<bool>("sphere")}, shapeFactor_{1},
      residualRe_{dict.lookupOrDefault<scalar>("residualRe", 1e-6)}
{
    if (residualRe_ <= 0) {
        FatalErrorInFunction << "Found invalid residualRe = " << residualRe_
                             << ", expected to be small positive number."
                             << exit(FatalError);
    }

    if (not sphere_) {
        shapeFactor_ = dict.lookup<scalar>("factor");
    }
}


Foam::tmp<Foam::volScalarField> Foam::heatTransferModels::Gnielinski::K(
    const scalar residualAlpha) const
{
    const volScalarField Re(max(interface_.Re(), residualRe_));
    const volScalarField Pr(interface_.Pr());

    auto Nu_laminar = 0.664 * sqrt(Re) * cbrt(Pr);
    auto Nu_turbulence = 0.037 * pow(Re, 0.8) * Pr /
        (1 + 2.443 * pow(Re, -0.1) * (pow(Pr, 2. / 3.) - 1));

    volScalarField Nu(IOobject("Nu_" + interface_.name(), interface_.mesh(),
                               IOobject::NO_READ, IOobject::AUTO_WRITE, true),
                      2 + sqrt(sqr(Nu_laminar) + sqr(Nu_turbulence)));

    if (interface_.mesh().time().writeTime()) {
        Nu.write();
    }
    
    if (sphere_) {
        const volScalarField &void_fraction = interface_.continuous();
        Nu *= (1 + 1.5 * (1 - void_fraction));
    }
    else {
        Nu *= shapeFactor_;
    }

    return 6 * max(interface_.dispersed(), residualAlpha) *
        interface_.continuous().thermo().kappa() * Nu /
        sqr(interface_.dispersed().d());
}