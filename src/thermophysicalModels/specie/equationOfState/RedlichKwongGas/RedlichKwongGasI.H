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

#include "RedlichKwongGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
inline Foam::RedlichKwongGas<Specie>::RedlichKwongGas(
    const Specie &sp, const realFluidProperty &prop)
    : Specie(sp), property_(prop)
{

    static_assert(Foam::is_cubic_EOS<Foam::RedlichKwongGas<Specie>>::value,
                  "RK must be cubic EOS");
}

template <class Specie>
inline Foam::RedlichKwongGas<Specie>::RedlichKwongGas(const word &name,
                                                      const RedlichKwongGas &pg)
    : Specie(name, pg), property_(pg.property_)
{}

template <class Specie>
inline Foam::autoPtr<Foam::RedlichKwongGas<Specie>> Foam::RedlichKwongGas<
    Specie>::clone() const
{
    return autoPtr<RedlichKwongGas<Specie>>(new RedlichKwongGas<Specie>(*this));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Specie>
inline void Foam::RedlichKwongGas<Specie>::operator+=(
    const RedlichKwongGas<Specie> &pg)
{
    noCoefficientMixing(RedlichKwongGas);
}

template <class Specie>
inline void Foam::RedlichKwongGas<Specie>::operator*=(const scalar s)
{
    noCoefficientMixing(RedlichKwongGas);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class Specie>
Foam::RedlichKwongGas<Specie> Foam::operator+(
    const RedlichKwongGas<Specie> &pg1, const RedlichKwongGas<Specie> &pg2)
{
    noCoefficientMixing(RedlichKwongGas);
}

template <class Specie>
Foam::RedlichKwongGas<Specie> Foam::operator*(const scalar s,
                                              const RedlichKwongGas<Specie> &pg)
{
    noCoefficientMixing(RedlichKwongGas);
}

template <class Specie>
Foam::RedlichKwongGas<Specie> Foam::operator==(
    const RedlichKwongGas<Specie> &pg1, const RedlichKwongGas<Specie> &pg2)
{
    noCoefficientMixing(RedlichKwongGas);
}
