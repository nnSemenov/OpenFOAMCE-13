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

#include "PengRobinson.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson(
    const Specie& sp,
    const realFluidProperty& prop
)
:
    Specie(sp),
    property_(prop)
{

}


template<class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson(
    const word& name,
                                                const PengRobinson &pg)
:
    Specie(name, pg),
    property_(pg.property_)
{
}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinson<Specie>> Foam::PengRobinson<
    Specie>::clone() const
{
    return autoPtr<PengRobinson<Specie>>(new PengRobinson<Specie>(*this));
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PengRobinson<Specie>::operator+=(
    const PengRobinson<Specie> &pg)
{
    noCoefficientMixing(PengRobinson);
}


template<class Specie>
inline void Foam::PengRobinson<Specie>::operator*=(const scalar s)
{
    noCoefficientMixing(PengRobinson);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::PengRobinson<Specie> Foam::operator+(const PengRobinson<Specie> &pg1,
                                           const PengRobinson<Specie> &pg2)
{
    noCoefficientMixing(PengRobinson);
}


template<class Specie>
Foam::PengRobinson<Specie> Foam::operator*(
    const scalar s,
                                           const PengRobinson<Specie> &pg)
{
    noCoefficientMixing(PengRobinson);
}


template<class Specie>
Foam::PengRobinson<Specie> Foam::operator==(const PengRobinson<Specie> &pg1,
                                            const PengRobinson<Specie> &pg2)
{
    noCoefficientMixing(PengRobinson);
}


// ************************************************************************* //
