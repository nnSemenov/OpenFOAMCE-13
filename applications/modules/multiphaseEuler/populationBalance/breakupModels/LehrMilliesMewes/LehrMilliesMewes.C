/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LehrMilliesMewes.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{
    defineTypeNameAndDebug(LehrMilliesMewes, 0);
    addToRunTimeSelectionTable(breakupModel, LehrMilliesMewes, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::breakupModels::LehrMilliesMewes::LehrMilliesMewes
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binary(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volInternalScalarField>
Foam::populationBalance::breakupModels::LehrMilliesMewes::rate
(
    const label i,
    const label j
) const
{
    using Foam::constant::mathematical::pi;

    const dimensionedScalar& dSphi = popBal_.dSph(i);
    const dimensionedScalar& dSphj = popBal_.dSph(j);

    const volInternalScalarField& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(j));
    const volInternalScalarField& sigma = tsigma();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volInternalScalarField& epsilonc = tepsilonc();

    volInternalScalarField L
    (
        pow(sigma/rhoc, static_cast<scalar>(3.0/5.0))/pow(epsilonc, static_cast<scalar>(2.0/5.0))
    );

    // Reset of dimension to pure length to avoid problems in transcendental
    // functions due to small exponents
    L.dimensions().reset(dimensions::length);

    const volInternalScalarField T
    (
        pow(sigma/rhoc, static_cast<scalar>(2.0/5.0))/pow(epsilonc, static_cast<scalar>(3.0/5.0))
    );

    return
        scalar{0.5}*pow(dSphj/L, static_cast<scalar>(5.0/3.0))
       *exp(-sqrt(scalar{2.0})/pow3(dSphj/L))
       *6/pow(pi, scalar{1.5})/pow3(dSphi/L)
       *exp(-static_cast<scalar>(9.0/4.0)*sqr(log(pow(scalar{2}, scalar{0.4})*dSphi/L)))
       /max(1 + erf(scalar{1.5}*log(pow(scalar{2.0}, static_cast<scalar>(1.0/15.0))*dSphj/L)), small)
       /(T*pow3(L));
}


// ************************************************************************* //
