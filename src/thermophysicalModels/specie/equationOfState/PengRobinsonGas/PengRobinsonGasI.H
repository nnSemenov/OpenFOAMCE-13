/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "PengRobinsonGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const Specie& sp,
    const realFluidProperty& prop
)
:
    Specie(sp),
    property_(prop)
{

}


template<class Specie>
inline Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const word& name,
    const PengRobinsonGas& pg
)
:
    Specie(name, pg),
    property_(pg.property_)
{
}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinsonGas <Specie>>
Foam::PengRobinsonGas<Specie>::clone() const
{
    return autoPtr<PengRobinsonGas<Specie>>
    (
        new PengRobinsonGas<Specie>(*this)
    );
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PengRobinsonGas<Specie>::operator+=
(
    const PengRobinsonGas<Specie>& pg
)
{
    noCoefficientMixing(RedlichKwongGas);
}


template<class Specie>
inline void Foam::PengRobinsonGas<Specie>::operator*=(const scalar s)
{
    noCoefficientMixing(RedlichKwongGas);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator+
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    noCoefficientMixing(RedlichKwongGas);
}


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator*
(
    const scalar s,
    const PengRobinsonGas<Specie>& pg
)
{
    noCoefficientMixing(RedlichKwongGas);
}


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator==
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    noCoefficientMixing(RedlichKwongGas);
}


// ************************************************************************* //
