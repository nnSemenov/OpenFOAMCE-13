/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    const scalar& Tc,
    const scalar& Vc,
    const scalar& Zc,
    const scalar& Pc,
    const scalar& omega
)
:
    Specie(sp),
    Tc_(Tc),
    Vc_(Vc),
    Zc_(Zc),
    Pc_(Pc),
    omega_(omega)
{}


template<class Specie>
inline Foam::cubicEOSCoefficient Foam::PengRobinsonGas<Specie>::coeff_AB(scalar p, scalar T) const {

    const scalar Tr = T/Tc_;
    const scalar Pr = p/Pc_;
    const scalar a = 0.45724*sqr(RR*Tc_)/Pc_;
    const scalar b = 0.07780*RR*Tc_/Pc_;
    const scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
    const scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));

    const scalar A = a*alpha*p/sqr(RR*T);
    const scalar B = b*p/(RR*T);
    return Foam::cubicEOSCoefficient{a,b,A,B,Tr,Pr,kappa,alpha};
}

template<class Specie>
inline Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const word& name,
    const PengRobinsonGas& pg
)
:
    Specie(name, pg),
    Tc_(pg.Tc_),
    Vc_(pg.Vc_),
    Zc_(pg.Zc_),
    Pc_(pg.Pc_),
    omega_(pg.omega_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinsonGas <Specie>>
Foam::PengRobinsonGas<Specie>::clone() const
{
    return autoPtr<PengRobinsonGas<Specie>>
    (
        new PengRobinsonGas<Specie>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return p/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::h(scalar p, scalar T) const
{
    const auto coeff = this->coeff_AB(p,T);
//    const scalar A=coeff.A;
    const scalar B=coeff.B;
    const scalar Tr=coeff.Tr;
    const scalar kappa = coeff.kappa;
    const scalar alpha = coeff.alpha;

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *Tc_
       *(
           Tr*(Z - 1)
         - 2.078*(1 + kappa)*sqrt(alpha)
          *log((Z + 2.414*B)/(Z - 0.414*B))
        );
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::Cp(scalar p, scalar T) const
{
    const auto coeff = this->coeff_AB(p,T);
    const scalar a=coeff.a;
    const scalar b=coeff.b;
    const scalar A=coeff.A;
    const scalar B=coeff.B;

    const scalar Z = this->Z(p, T);

    const scalar kappa = coeff.kappa;
    const scalar ap = kappa*a*(kappa/Tc_ - (1 + kappa)/sqrt(T*Tc_));
    const scalar app = kappa*a*(1 + kappa)/(2*sqrt(pow3(T)*Tc_));

    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b*RR);

    const scalar root2 = sqrt(2.0);

    return
    (
        app*(T/(2*root2*b))*log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
      + RR*sqr(M - N)/(sqr(M) - 2*A*(Z + B))
      - RR
    )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::e(scalar p, scalar T) const
{
    const auto coeff = this->coeff_AB(p,T);
    const scalar B=coeff.B;
    const scalar kappa = coeff.kappa;
    const scalar alpha = coeff.alpha;

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *Tc_
       *(
         - 2.078*(1 + kappa)*sqrt(alpha)
          *log((Z + 2.414*B)/(Z - 0.414*B))
        );
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::Cv(scalar p, scalar T) const
{
    const auto coeff = this->coeff_AB(p,T);
    const scalar a=coeff.a;
    const scalar b=coeff.b;
    const scalar B=coeff.B;
    const scalar kappa = coeff.kappa;
    const scalar Z = this->Z(p, T);

    const scalar app = kappa*a*(1 + kappa)/(2*sqrt(pow3(T)*Tc_));

    const scalar root2 = sqrt(2.0);

    return
    (
        app*(T/(2*root2*b))*log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
      - RR
    )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::sp
(
    scalar p,
    scalar T
) const
{
    const auto coeff = this->coeff_AB(p,T);
    const scalar Tr = coeff.Tr;
    const scalar B = coeff.B;
    const scalar kappa = coeff.kappa;

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *(
          - log(p/Pstd)
          + (
                log(Z - B)
              - 2.078*kappa*((1 + kappa)/sqrt(Tr) - kappa)
               *log((Z + 2.414*B)/(Z - 0.414*B))
            )
        );
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::sv
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);

    return 1.0/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const auto coeff = this->coeff_AB(p,T);
    const scalar A=coeff.A;
    const scalar B=coeff.B;

    const scalar a2 = B - 1;
    const scalar a1 = A - 2*B - 3*sqr(B);
    const scalar a0 = -A*B + sqr(B) + pow3(B);

    return solveCubicEquation(a2, a1, a0);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    const auto coeff = this->coeff_AB(p,T);
    const scalar A=coeff.A;
    const scalar B=coeff.B;
    const scalar kappa = coeff.kappa;
    const scalar a=coeff.a;
    const scalar b=coeff.b;

    const scalar Z = this->Z(p, T);

    const scalar ap = kappa*a*(kappa/Tc_ - (1 + kappa)/sqrt(T*Tc_));
    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b*RR);

    return this->R()*sqr(M - N)/(sqr(M) - 2*A*(Z + B));
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::alphav
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PengRobinsonGas<Specie>::operator+=
(
    const PengRobinsonGas<Specie>& pg
)
{
    scalar X1 = this->Y()/this->W();
    Specie::operator+=(pg);

    if (mag(this->Y()) > small)
    {
        X1 *= this->W()/this->Y();
        const scalar X2 = this->W()*pg.Y()/(pg.W()*this->Y());

        Tc_ = X1*Tc_ + X2*pg.Tc_;
        Vc_ = X1*Vc_ + X2*pg.Vc_;
        Zc_ = X1*Zc_ + X2*pg.Zc_;
        Pc_ = RR*Zc_*Tc_/Vc_;
        omega_ = X1*omega_ + X2*pg.omega_;
    }
}


template<class Specie>
inline void Foam::PengRobinsonGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator+
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return PengRobinsonGas<Specie>
        (
            sp,
            pg1.Tc_,
            pg1.Vc_,
            pg1.Zc_,
            pg1.Pc_,
            pg1.omega_
        );
    }
    else
    {
        const scalar X1 = sp.W()*pg1.Y()/(pg1.W()*sp.Y());
        const scalar X2 = sp.W()*pg2.Y()/(pg2.W()*sp.Y());

        const scalar Tc = X1*pg1.Tc_ + X2*pg2.Tc_;
        const scalar Vc = X1*pg1.Vc_ + X2*pg2.Vc_;
        const scalar Zc = X1*pg1.Zc_ + X2*pg2.Zc_;

        return PengRobinsonGas<Specie>
        (
            sp,
            Tc,
            Vc,
            Zc,
            RR*Zc*Tc/Vc,
            X1*pg1.omega_ + X2*pg2.omega_
        );
    }
}


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator*
(
    const scalar s,
    const PengRobinsonGas<Specie>& pg
)
{
    return PengRobinsonGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Zc_,
        pg.Pc_,
        pg.omega_
    );
}


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator==
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    return PengRobinsonGas<Specie>
    (
        static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2),
        NaN,
        NaN,
        NaN,
        NaN,
        NaN
    );
}


// ************************************************************************* //
