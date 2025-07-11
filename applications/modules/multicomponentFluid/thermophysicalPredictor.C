/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "multicomponentFluid.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multicomponentFluid::thermophysicalPredictor()
{
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi_,
            mesh.schemes().div("div(phi,Yi_h)")
        )
    );

    reaction->correct();

    forAll(Y(), i)
    {
        volScalarField& Yi = Y_()[i];

        if (thermo_().solveSpecie(i))
        {
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho_, Yi)
              + mvConvection->fvmDiv(phi_, Yi)
              + thermophysicalTransport->divj(Yi)
             ==
                reaction->R(Yi)
              + fvModels().source(rho_, Yi)
            );

            YiEqn.relax();

            fvConstraints().constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints().constrain(Yi);
        }
        else
        {
            Yi.correctBoundaryConditions();
        }
    }

    thermo_().normaliseY();


    volScalarField& he = thermo_().he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho_, he) + mvConvection->fvmDiv(phi_, he)
      + fvc::ddt(rho_, K) + fvc::div(phi_, K)
      + pressureWork
        (
            he.name() == "e"
          ? mvConvection->fvcDiv(phi_, p()/rho_)()
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        reaction->Qdot()
      + (
            buoyancy.valid()
          ? fvModels().source(rho_, he) + rho_*(U_ & buoyancy->g)
          : fvModels().source(rho_, he)
        )
    );

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(he);

    thermo_().correct();
}


// ************************************************************************* //
