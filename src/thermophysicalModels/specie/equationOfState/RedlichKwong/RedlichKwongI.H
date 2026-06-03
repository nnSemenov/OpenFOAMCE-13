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

#include "RedlichKwong.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
inline Foam::RedlichKwong<Specie>::RedlichKwong(const Specie &sp, const realFluidProperty &prop)
    : Specie(sp), property_(prop)
{

    static_assert(Foam::is_cubic_EOS<Foam::RedlichKwong<Specie>>::value,
                  "RK must be cubic EOS");
}

template <class Specie>
inline Foam::RedlichKwong<Specie>::RedlichKwong(const word &name,
                                                const RedlichKwong &pg)
    : Specie(name, pg), property_(pg.property_)
{}

template <class Specie>
inline Foam::autoPtr<Foam::RedlichKwong<Specie>> Foam::RedlichKwong<
    Specie>::clone() const
{
    return autoPtr<RedlichKwong<Specie>>(new RedlichKwong<Specie>(*this));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Specie>
inline void Foam::RedlichKwong<Specie>::operator+=(
    const RedlichKwong<Specie> &pg)
{
    noCoefficientMixing(RedlichKwong);
}

template <class Specie>
inline void Foam::RedlichKwong<Specie>::operator*=(const scalar s)
{
    noCoefficientMixing(RedlichKwong);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class Specie>
Foam::RedlichKwong<Specie> Foam::operator+(const RedlichKwong<Specie> &pg1,
                                           const RedlichKwong<Specie> &pg2)
{
    noCoefficientMixing(RedlichKwong);
}

template <class Specie>
Foam::RedlichKwong<Specie> Foam::operator*(const scalar s,
                                           const RedlichKwong<Specie> &pg)
{
    noCoefficientMixing(RedlichKwong);
}

template <class Specie>
Foam::RedlichKwong<Specie> Foam::operator==(const RedlichKwong<Specie> &pg1,
                                            const RedlichKwong<Specie> &pg2)
{
    noCoefficientMixing(RedlichKwong);
}
