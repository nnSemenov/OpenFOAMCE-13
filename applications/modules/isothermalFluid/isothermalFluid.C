/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "isothermalFluid.H"
#include "localEulerDdtScheme.H"
#include "hydrostaticInitialisation.H"
#include "fvcMeshPhi.H"
#include "fvcVolumeIntegrate.H"
#include "fvcReconstruct.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(isothermalFluid, 0);
    addToRunTimeSelectionTable(solver, isothermalFluid, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::isothermalFluid::correctCoNum()
{
    fluidSolver::correctCoNum(rho_, phi_);
}


void Foam::solvers::isothermalFluid::continuityErrors()
{
    fluidSolver::continuityErrors(rho_, thermo().rho(), phi_);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::solvers::isothermalFluid::pressureWork
(
    const tmp<volScalarField::Internal>& work
) const
{
    if (mesh.moving())
    {
        return
            work
          + fvc::div
            (
                fvc::interpolate(rho_)*fvc::meshPhi(rho_, U_),
                p()/rho_,
                "div(phi,(p|rho))"
            )();
    }
    else
    {
        return move(work);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::isothermalFluid::isothermalFluid
(
    fvMesh& mesh,
    autoPtr<fluidThermo> thermoPtr
)
:
    fluidSolver(mesh),

    thermoPtr_(thermoPtr),

    rho_
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermoPtr_().renameRho()
    ),

    dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.name(),
            mesh
        ),
        mesh,
        dimensionedScalar(p_().dimensions()/dimTime, 0)
    ),

    buoyancy(buoyancy::New(mesh)),

    pressureReference
    (
        p_(),
        p_rgh_(),
        pimple.dict(),
        thermoPtr_().incompressible()
    ),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho_*U_) & mesh.Sf()
    ),

    K("K", 0.5*magSqr(U_)),

    momentumTransport
    (
        compressible::momentumTransportModel::New
        (
            rho_,
            U_,
            phi_,
            thermoPtr_()
        )
    ),

    initialMass(fvc::domainIntegrate(rho_)),

    MRF(mesh)

//    thermo(thermoPtr_()),
//    p(p_),
//    p_rgh(p_rgh_),
//    rho(rho_),
//    U(U_),
//    phi(phi_)
{
    mesh.schemes().setFluxRequired(p_().name());
    momentumTransport->validate();

    if (buoyancy.valid())
    {
        hydrostaticInitialisation
        (
            p_rgh_(),
            p_(),
            rho_,
            U_,
            buoyancy->gh,
            buoyancy->ghf,
            buoyancy->pRef,
            thermoPtr_(),
            pimple.dict()
        );

        netForce = new volVectorField
        (
            IOobject
            (
                "netForce",
                runTime.name(),
                mesh
            ),
            fvc::reconstruct
            (
                (-buoyancy->ghf*fvc::snGrad(rho_) - fvc::snGrad(p_rgh_()))
               *mesh.magSf()
            )
        );
    }

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


Foam::solvers::isothermalFluid::isothermalFluid(fvMesh& mesh)
:
    isothermalFluid(mesh, fluidThermo::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::isothermalFluid::~isothermalFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::isothermalFluid::preSolve()
{
    if ((mesh.dynamic() || MRF.size()) && !rhoUf.valid())
    {
        Info<< "Constructing face momentum rhoUf" << endl;

        // Ensure the U BCs are up-to-date before constructing Uf
        U_.correctBoundaryConditions();

        rhoUf = new surfaceVectorField
        (
            IOobject
            (
                "rhoUf",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(rho_*U_)
        );
    }

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Store divrhoU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    if (correctPhi || mesh.topoChanging())
    {
        divrhoU = new volScalarField
        (
            "divrhoU",
            fvc::div(fvc::absolute(phi_, rho_, U_))
        );
    }

    fvModels().preUpdateMesh();

    // Store momentum to set rhoUf for introduced faces
    if (mesh.topoChanging())
    {
        rhoU = new volVectorField("rhoU", rho_*U_);

        for (label i = 1; i <= rhoUf().nOldTimes(false); ++ i)
        {
            rhoU().oldTimeRef(i) == rho_.oldTime(i)*U_.oldTime(i);
        }
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::isothermalFluid::prePredictor()
{
    if
    (
        !mesh.schemes().steady()
     && !pimple.simpleRho()
     && pimple.firstIter()
    )
    {
        correctDensity();
    }
}


void Foam::solvers::isothermalFluid::momentumTransportPredictor()
{
    momentumTransport->predict();
}


void Foam::solvers::isothermalFluid::thermophysicalTransportPredictor()
{}


void Foam::solvers::isothermalFluid::thermophysicalPredictor()
{
    thermoPtr_().correct();
}


void Foam::solvers::isothermalFluid::pressureCorrector()
{
    while (pimple.correct())
    {
        if (buoyancy.valid())
        {
            correctBuoyantPressure();
        }
        else
        {
            correctPressure();
        }
    }

    tUEqn.clear();
}


void Foam::solvers::isothermalFluid::momentumTransportCorrector()
{
    momentumTransport->correct();
}


void Foam::solvers::isothermalFluid::thermophysicalTransportCorrector()
{}


void Foam::solvers::isothermalFluid::postSolve()
{
    divrhoU.clear();

    if (!mesh.schemes().steady())
    {
        rho_ = thermo().rho();

        // Correct rhoUf with the updated density if the mesh is moving
        fvc::correctRhoUf(rhoUf, rho_, U_, phi_, MRF);
    }
}


// ************************************************************************* //
