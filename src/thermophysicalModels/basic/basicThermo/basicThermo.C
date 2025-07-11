/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "basicThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "gradientEnergyCalculatedTemperatureFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "mixedEnergyCalculatedTemperatureFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "energyJumpFvPatchScalarField.H"
#include "energyFvScalarFieldSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicThermo, 0);
    defineRunTimeSelectionTable(basicThermo, fvMesh);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicThermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const char* name
)
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
    {
        volScalarField* fPtr
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volScalarField>(name);
}


Foam::wordList Foam::basicThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");

            // If the number of number of components in the name
            // is greater than nCmpt return an empty list
            if (i == nCmpt)
            {
                return wordList();
            }
        }
        beg = end + 1;
    }

    // If the number of number of components in the name is not equal to nCmpt
    // return an empty list
    if (i + 1 != nCmpt)
    {
        return wordList();
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i].replaceAll(">","");
    }

    return cmpts;
}


Foam::List<Foam::Pair<Foam::word>> Foam::basicThermo::thermoNameComponents
(
    const word& thermoName
)
{
    const wordList components(splitThermoName(thermoName, 5));

    return List<Pair<word>>
    {
        {"transport", components[0]},
        {"thermo", components[1]},
        {"equationOfState", components[2]},
        {"specie", components[3]},
        {"energy", components[4]}
    };
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basicThermo::heBoundaryBaseTypes()
{
    const volScalarField::Boundary& tbf = T().boundaryField();

    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        if (tbf[patchi].overridesConstraint())
        {
            hbt[patchi] = tbf[patchi].patch().type();
        }
    }

    return hbt;
}


Foam::wordList Foam::basicThermo::heBoundaryTypes()
{
    const volScalarField::Boundary& tbf = T().boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
         || isA<gradientEnergyCalculatedTemperatureFvPatchScalarField>
            (
                tbf[patchi]
            )
        )
        {
            hbt[patchi] = gradientEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<mixedFvPatchScalarField>(tbf[patchi])
         || isA<mixedEnergyCalculatedTemperatureFvPatchScalarField>
            (
                tbf[patchi]
            )
        )
        {
            hbt[patchi] = mixedEnergyFvPatchScalarField::typeName;
        }
        else if (isA<jumpCyclicFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpFvPatchScalarField::typeName;
        }
    }

    return hbt;
}


Foam::HashTable<Foam::word> Foam::basicThermo::heSourcesTypes()
{
    const HashTable<word> tst = T().sources().types();

    HashTable<word> hst;
    forAllConstIter(typename HashTable<word>, tst, iter)
    {
        hst.set(iter.key(), energyFvScalarFieldSource::typeName);
    }

    return hst;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicThermo::implementation::implementation
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    mesh_(mesh),

    phaseName_(phaseName),

    T_
    (
        IOobject
        (
            phasePropertyName("T", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    kappa_
    (
        IOobject
        (
            phasePropertyName("kappa", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimThermalConductivity, Zero)
    ),

    dpdt_(dict.lookupOrDefault<Switch>("dpdt", true)),

    pOffset_(dict.lookupOrDefault<dimensionedScalar>(
            "pOffset",
            dimensionedScalar{"pOffset",dimPressure,0.0}
    ))
{

    Info<<"Pressure offset is "<<this->pOffset().value()<<" Pa"<<endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicThermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return New<basicThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermo::~basicThermo()
{}


Foam::basicThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicThermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!(he().name() == phasePropertyName(a)))
    {
        FatalErrorInFunction
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


Foam::tmp<Foam::volScalarField> Foam::basicThermo::gamma() const
{
    return volScalarField::New(phasePropertyName("gamma"), Cp()/Cv());
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp(T, patchi)/Cv(T, patchi);
}


const Foam::volScalarField& Foam::basicThermo::implementation::T() const
{
    return T_;
}


Foam::volScalarField& Foam::basicThermo::implementation::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicThermo::implementation::kappa() const
{
    return kappa_;
}


void Foam::basicThermo::implementation::read(const dictionary&dict)
{}


const Foam::dimensionedScalar & Foam::basicThermo::implementation::pOffset() const {
    return this->pOffset_;
}

Foam::dimensionedScalar & Foam::basicThermo::implementation::pOffset() {
    return this->pOffset_;
}

// ************************************************************************* //
