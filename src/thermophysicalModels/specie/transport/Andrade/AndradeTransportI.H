/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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

#include "AndradeTransport.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::AndradeTransport<Thermo>::AndradeTransport
(
    const Thermo& t,
    const coeffList& muCoeffs,
    const coeffList& kappaCoeffs
)
:
    Thermo(t),
    muCoeffs_(muCoeffs),
    kappaCoeffs_(kappaCoeffs)
{}


template<class Thermo>
inline Foam::AndradeTransport<Thermo>::AndradeTransport
(
    const word& name,
    const AndradeTransport& pt
)
:
    Thermo(name, pt),
    muCoeffs_(pt.muCoeffs_),
    kappaCoeffs_(pt.kappaCoeffs_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::AndradeTransport<Thermo>>
Foam::AndradeTransport<Thermo>::clone() const
{
    return autoPtr<AndradeTransport<Thermo>>
    (
        new AndradeTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::AndradeTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return exp
    (
        muCoeffs_[0]
      + muCoeffs_[1]*T
      + muCoeffs_[2]*sqr(T)
      + muCoeffs_[3]/(muCoeffs_[4] + T)
    );
}


template<class Thermo>
inline Foam::scalar Foam::AndradeTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return exp
    (
        kappaCoeffs_[0]
      + kappaCoeffs_[1]*T
      + kappaCoeffs_[2]*sqr(T)
      + kappaCoeffs_[3]/(kappaCoeffs_[4] + T)
    );
}

template <class Thermo>
void Foam::AndradeTransport<Thermo>::operator+=(const AndradeTransport<Thermo>& at) {

    scalar Y1 = this->Y();

    Thermo::operator+=(at);
    // this->Y() gives changed results after this +=
    if (mag(this->Y())>small) {
        Y1/=this->Y();
        const scalar Y2=at.Y()/this->Y();
        forAll(muCoeffs_,idx) {
            muCoeffs_[idx]=Y1*this->muCoeffs_[idx] + Y2*at.muCoeffs_[idx];
        }
        forAll(kappaCoeffs_, idx) {
            kappaCoeffs_[idx]=Y1*this->kappaCoeffs_[idx] + Y2*at.kappaCoeffs_[idx];
        }
    }
}

template <class Thermo>
void Foam::AndradeTransport<Thermo>::operator*=(const scalar s) {
    Thermo::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::AndradeTransport<Thermo> Foam::operator+(
    const AndradeTransport<Thermo>&at1,
    const AndradeTransport<Thermo>& at2
    ) {
    Thermo t(
        static_cast<const Thermo&>(at1) + static_cast<const Thermo&>(at2));

    if (mag(t.Y())<small) {
        return AndradeTransport<Thermo> (
            t,
            0,
            at1.muCoeffs_,
            at1.kappaCoeffs_
        );
    }
    const scalar Y1=at1.Y()/t.Y();
    const scalar Y2=at2.Y()/t.Y();
    decltype(at1.muCoeffs_) muCoeff, kappaCoeff;
    forAll(muCoeff, idx) {
        muCoeff[idx]    = Y1*at1.muCoeffs_[idx]      +Y2*at2.muCoeffs_[idx];
    }
    forAll(kappaCoeff, idx) {
        kappaCoeff[idx] = Y1*at1.kappaCoeffs_[idx]   +Y2*at2.kappaCoeffs_[idx];
    }
    return AndradeTransport<Thermo> (
        t,
        muCoeff,
        kappaCoeff
    );
}

template<class Thermo>
inline Foam::AndradeTransport<Thermo> Foam::operator*(
    const scalar s,
    const AndradeTransport<Thermo>& at
) {
    return AndradeTransport<Thermo> (
        s*static_cast<const Thermo&>(at),
        at.muCoeffs_,
        at.kappaCoeffs_
    );
}

// ************************************************************************* //
