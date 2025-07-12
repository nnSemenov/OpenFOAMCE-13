/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "fluid.H"
#include "fvcDdt.H"
#include "fvmDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::fluid::thermophysicalPredictor()
{
    volScalarField& he = this->thermoPtr_().he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho_, he) + fvm::div(phi_, he)
      + fvc::ddt(rho_, K) + fvc::div(phi_, K)
      + pressureWork
        (
            he.name() == "e"
          ? fvc::div(phi_, p_()/rho_)()
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        (
            buoyancy.valid()
          ? fvModels().source(rho_, he) + rho_*(U_ & buoyancy->g)
          : fvModels().source(rho_, he)
        )
    );

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(he);

    this->thermoPtr_().correct();
}


// ************************************************************************* //
