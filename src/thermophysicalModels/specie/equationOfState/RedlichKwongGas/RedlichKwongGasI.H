
#include "RedlichKwongGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::RedlichKwongGas<Specie>::RedlichKwongGas
        (
                const Specie &sp,
                const scalar &Tc,
                const scalar &Vc,
                const scalar &Zc,
                const scalar &Pc
        )
        :
        Specie(sp),
        Tc_(Tc),
        Vc_(Vc),
        Zc_(Zc),
        Pc_(Pc) {}


template<class Specie>
inline Foam::cubicEOSCoefficient Foam::RedlichKwongGas<Specie>::coeff_AB(scalar p, scalar T) const {

    const scalar Tr = T / Tc_;
    const scalar Pr = p / Pc_;
    const scalar a = 0.42748 * sqr(RR * Tc_) * sqrt(Tc_) / Pc_;
    const scalar b = 0.08664 * RR * Tc_ / Pc_;

    const scalar A = a / p / sqrt(T) * sqr(p / (RR * T));
    const scalar B = b * p / (RR * T);
    return Foam::cubicEOSCoefficient{a, b, A, B, Tr, Pr, 0, 1};
}

template<class Specie>
inline Foam::RedlichKwongGas<Specie>::RedlichKwongGas
        (
                const word &name,
                const RedlichKwongGas &pg
        )
        :
        Specie(name, pg),
        Tc_(pg.Tc_),
        Vc_(pg.Vc_),
        Zc_(pg.Zc_),
        Pc_(pg.Pc_) {}


template<class Specie>
inline Foam::autoPtr<Foam::RedlichKwongGas<Specie>>
Foam::RedlichKwongGas<Specie>::clone() const {
    return autoPtr<RedlichKwongGas<Specie>>
            (
                    new RedlichKwongGas<Specie>(*this)
            );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::rho
        (
                scalar p,
                scalar T
        ) const {
    const scalar Z = this->Z(p, T);
    return p / (Z * this->R() * T);
}

//template<class Specie>
//inline Foam::scalar Foam::RedlichKwongGas<Specie>::partialA_partialT_const_p(const scalar p, const scalar T,
//                                                                             Foam::cubicEOSCoefficient &coeff) const {
//    const scalar a = coeff.a;
//    return -2.5 * a * p / sqr(RR) / pow(T, 3.5);
//}
//
//template<class Specie>
//inline Foam::scalar
//Foam::RedlichKwongGas<Specie>::partialZ_partialT_const_p(const scalar p, const scalar T, const scalar Z,
//                                                         Foam::cubicEOSCoefficient &coeff) const {
//    const auto parAparT_p = this->partialA_partialT_const_p(p, T, coeff);
//    const scalar A = coeff.A;
//    const scalar B = coeff.B;
//    assert(Z > 0);
//
//    return parAparT_p * (B - Z) / (3 * sqr(Z) - 2 * Z + (A - B - sqr(B)));
//}

template<class Specie>
inline Foam::residualEnergy Foam::RedlichKwongGas<Specie>::residualHE(scalar p, scalar T) const {

    const auto coeff = this->coeff_AB(p, T);
    const scalar A = coeff.A;
    const scalar B = coeff.B;
    const scalar Z = this->Z(p, T);

    const scalar e_R_dimless = - (3 * A) / (2 * B) * log(Z/(Z+B));
    const scalar h_R_dimless =
            (Z - 1)
            + e_R_dimless;
    const scalar e_R = T * this->R() * e_R_dimless;
    const scalar h_R = T * this->R() * h_R_dimless;
    return Foam::residualEnergy{h_R, e_R, Z};
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::h(scalar p, scalar T) const {
//    return 0;
    return this->residualHE(p, T).h_R;
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::e(scalar p, scalar T) const {
//    return 0;
    return this->residualHE(p, T).e_R;
}

template<class Specie>
inline Foam::residualSpecificHeat
Foam::RedlichKwongGas<Specie>::residualSpecificHeat(const scalar p, const scalar T) const {

    const auto coeff = this->coeff_AB(p, T);
    const scalar A = coeff.A;
    const scalar B = coeff.B;
    const scalar Z = this->Z(p, T);

    const scalar Cv_dimless = (3 * A) / (4 * B) * log(1+B/Z);

    const scalar CpRMCvR_dimless =
            -1 -
            sqr(1 / (Z - B) + A / (2 * Z * (Z + B)))
            /
            (-1 / sqr(Z - B) + A * (2 * Z + B) / sqr(Z) / sqr(Z + B));
    const scalar Cp_dimless = Cv_dimless + CpRMCvR_dimless;

    return Foam::residualSpecificHeat{
            Cp_dimless * this->R(),
            Cv_dimless * this->R(),
            CpRMCvR_dimless * this->R(),
            Z
    };
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::Cp(scalar p, scalar T) const {
//    return 0;
    return this->residualSpecificHeat(p, T).Cp_R;
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::Cv(scalar p, scalar T) const {
//    return 0;
    return this->residualSpecificHeat(p, T).Cv_R;
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::CpMCv
        (
                scalar p,
                scalar T
        ) const {
    return this->residualSpecificHeat(p, T).Cp_R_minus_Cv_R;
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::sp
        (
                scalar p,
                scalar T
        ) const {
    NotImplemented;
    const auto coeff = this->coeff_AB(p, T);
    const scalar Pr = coeff.Pr;
    const scalar Tr = coeff.Tr;
    const scalar B = coeff.B;
    const scalar kappa = coeff.kappa;
    const scalar alpha = coeff.alpha;

    const scalar Z = this->Z(p, T);

    return
            this->R()
            * (
                    -log(p / Pstd)
                    + (
                            log(Z - B)
                            - 2.078 * kappa * ((1 + kappa) / sqrt(Tr) - kappa)
                              * log((Z + 2.414 * B) / (Z - 0.414 * B))
                    )
            );
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::sv
        (
                scalar p,
                scalar T
        ) const {
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::psi
        (
                scalar p,
                scalar T
        ) const {
    const scalar Z = this->Z(p, T);

    return 1.0 / (Z * this->R() * T);
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::Z
        (
                scalar p,
                scalar T
        ) const {
    const auto coeff = this->coeff_AB(p, T);
    const scalar A = coeff.A;
    const scalar B = coeff.B;

    const scalar a2 = -1;
    const scalar a1 = A - B - sqr(B);
    const scalar a0 = -A * B;

    return solveCubicEquation(a2, a1, a0);
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongGas<Specie>::alphav
        (
                scalar p,
                scalar T
        ) const {
    NotImplemented;
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::RedlichKwongGas<Specie>::operator+=
        (
                const RedlichKwongGas<Specie> &pg
        ) {
    scalar X1 = this->Y() / this->W();
    Specie::operator+=(pg);

    if (mag(this->Y()) > small) {
        X1 *= this->W() / this->Y();
        const scalar X2 = this->W() * pg.Y() / (pg.W() * this->Y());

        Tc_ = X1 * Tc_ + X2 * pg.Tc_;
        Vc_ = X1 * Vc_ + X2 * pg.Vc_;
        Zc_ = X1 * Zc_ + X2 * pg.Zc_;
        Pc_ = RR * Zc_ * Tc_ / Vc_;
//        omega_ = X1*omega_ + X2*pg.omega_;
    }
}


template<class Specie>
inline void Foam::RedlichKwongGas<Specie>::operator*=(const scalar s) {
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::RedlichKwongGas<Specie> Foam::operator+
        (
                const RedlichKwongGas<Specie> &pg1,
                const RedlichKwongGas<Specie> &pg2
        ) {
    Specie sp
            (
                    static_cast<const Specie &>(pg1)
                    + static_cast<const Specie &>(pg2)
            );

    if (mag(sp.Y()) < small) {
        return RedlichKwongGas<Specie>
                (
                        sp,
                        pg1.Tc_,
                        pg1.Vc_,
                        pg1.Zc_,
                        pg1.Pc_,
                        pg1.omega_
                );
    } else {
        const scalar X1 = sp.W() * pg1.Y() / (pg1.W() * sp.Y());
        const scalar X2 = sp.W() * pg2.Y() / (pg2.W() * sp.Y());

        const scalar Tc = X1 * pg1.Tc_ + X2 * pg2.Tc_;
        const scalar Vc = X1 * pg1.Vc_ + X2 * pg2.Vc_;
        const scalar Zc = X1 * pg1.Zc_ + X2 * pg2.Zc_;

        return RedlichKwongGas<Specie>
                (
                        sp,
                        Tc,
                        Vc,
                        Zc,
                        RR * Zc * Tc / Vc
                );
    }
}


template<class Specie>
Foam::RedlichKwongGas<Specie> Foam::operator*
        (
                const scalar s,
                const RedlichKwongGas<Specie> &pg
        ) {
    return RedlichKwongGas<Specie>
            (
                    s * static_cast<const Specie &>(pg),
                    pg.Tc_,
                    pg.Vc_,
                    pg.Zc_,
                    pg.Pc_
            );
}


template<class Specie>
Foam::RedlichKwongGas<Specie> Foam::operator==
        (
                const RedlichKwongGas<Specie> &pg1,
                const RedlichKwongGas<Specie> &pg2
        ) {
    return RedlichKwongGas<Specie>
            (
                    static_cast<const Specie &>(pg1) == static_cast<const Specie &>(pg2),
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN
            );
}

