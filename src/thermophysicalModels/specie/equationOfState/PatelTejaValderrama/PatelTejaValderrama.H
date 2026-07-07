/*---------------------------------------------------------------------------*\
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

Class
    Foam::PatelTejaValderrama

Description
    PatelTejaValderrama cubic equation of state for gases and liquids. Mixing
    rule is missing currently. Will enable for realGasMulticomponentMixture
    later.

Usage
    \table
        Property     | Description
        Tc           | Critical temperature [K]
        Vc           | Critical volume [m^3/kmol]
        Pc           | Critical pressure [Pa]
        omega        | Acentric factor [-]
        phase        | Preferred phase (vapor/liquid)
    \endtable

    Example specification of the PatelTejaValderrama equation of state for
    methane:
    \verbatim
    equationOfState
    {
        Tc        190.55;
        Vc        0.0285;
        Pc        4.595e6;
        omega     0.0115;
        phase     vapor;
    }
    \endverbatim

SourceFiles
    PatelTejaValderrama.C

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_PATELTEJAVALDERRAMA_H
#define MIKENO_PATELTEJAVALDERRAMA_H

#include "specie.H"
#include "autoPtr.H"
#include "speciesTable.H"
#include "mathematicalConstants.H"
#include "../equationOfState/cubicEOS.H"
#include <cassert>
#include <cmath>

namespace Foam
{

// Forward declaration of friend functions and operators

template <class Specie>
class PatelTejaValderrama;

template <class Specie>
Ostream &operator<<(Ostream &, const PatelTejaValderrama<Specie> &);

template <class Specie>
PatelTejaValderrama<Specie> operator+(const PatelTejaValderrama<Specie> &,
                                      const PatelTejaValderrama<Specie> &);

template <class Specie>
PatelTejaValderrama<Specie> operator*(const scalar s,
                                      const PatelTejaValderrama<Specie> &);

template <class Specie>
PatelTejaValderrama<Specie> operator==(const PatelTejaValderrama<Specie> &,
                                       const PatelTejaValderrama<Specie> &);

struct PTVDimlessCoefficient
{
    scalar A;
    scalar B;
    scalar C;
};

template <class Specie>
class PatelTejaValderrama : public Specie
{
  protected:
    realFluidProperty property_;

  public:
    struct EOSCore
    {
        scalar a;
        scalar b;
        scalar c;
        scalar da_dT;
        scalar dd2a_dT2;
        bool is_gas;

        [[nodiscard]] inline PTVDimlessCoefficient dimlessCoeff(scalar p,
                                                                scalar T) const
        {
            const scalar A = a * p / sqr(RR * T);
            const scalar B = b * p / (RR * T);
            const scalar C = c * p / (RR * T);
            return {A, B, C};
        }

        [[nodiscard]] inline scalar R(scalar W) const { return RR / W; }

        [[nodiscard]] inline scalar Z(const PTVDimlessCoefficient &coeff) const
        {
            const auto [A, B, C] = coeff;
            const scalar a2 = C - 1;
            const scalar a1 = A - 2 * B * C - sqr(B) - (B + C);
            const scalar a0 = sqr(B) * C + B * C - A * B;
            assert(B > 0);
            const scalar Z = solveCubicEquation(a2, a1, a0, this->is_gas, B);
            assert(Z > B);
            return Z;
        }

        [[nodiscard]] inline scalar Z(scalar p, scalar T) const
        {
            auto coeff = this->dimlessCoeff(p, T);
            return this->Z(coeff);
        }
        [[nodiscard]] inline scalar rho(scalar p, scalar T, scalar W) const
        {
            const scalar Z = this->Z(p, T);
            return p / (Z * this->R(W) * T);
        }

        [[nodiscard]] inline Foam::residualEnergy residualHE(scalar p, scalar T,
                                                             scalar W) const
        {
            const scalar z = this->Z(p, T);
            const scalar V = z * RR * T / p;
            const scalar cbcb = sqrt(sqr(c) + 6 * b * c + sqr(b));
            const scalar aT_part = (da_dT * T - a) / (RR * T * cbcb);
            const scalar log_part =
                log((2 * V - cbcb + c + b) / (2 * V + cbcb + c + b));
            const scalar dimless_U_R = aT_part * log_part;
            const scalar dimless_H_R = dimless_U_R + z - 1;

            return {
                .h_R = dimless_H_R * T * this->R(W),
                .e_R = dimless_U_R * T * this->R(W),
                .Z = z,
            };
        }

        [[nodiscard]] inline Foam::residualSpecificHeat residualSpecificHeat(
            scalar p, scalar T, scalar W) const
        {
            const scalar z = this->Z(p, T);
            const scalar V = z * RR * T / p;
            const scalar cbcb = sqrt(sqr(c) + 6 * b * c + sqr(b));

            const scalar aT_part = -dd2a_dT2 * T / (RR * cbcb);
            const scalar log_part =
                log((2 * V - cbcb + c + b) / (2 * V + cbcb + c + b));

            const scalar dimless_CvR = aT_part * log_part;

            const scalar dimless_CpMCv_R =
                -1 + T / RR * dpdT_V(p, T, z) * dVdT_p(p, T, z);

            const scalar dimless_CpR = dimless_CpMCv_R + dimless_CvR;

            return {
                .Cp_R = dimless_CpR * this->R(W),
                .Cv_R = dimless_CvR * this->R(W),
                .Cp_R_minus_Cv_R = dimless_CpMCv_R * this->R(W),
                .Z = z,
            };
        }

        [[nodiscard]] inline Foam::residualSpecificHeatDerivative
        residualSpecificHeatDerivative(scalar p, scalar T, scalar W) const
        {
            const scalar z = this->Z(p, T);
            const scalar V = z * RR * T / p;
            const scalar cbcb = sqrt(sqr(c) + 6 * b * c + sqr(b));

            const scalar a_log_part = -dd2a_dT2 / cbcb *
                log((2 * V - cbcb + c + b) / (2 * V + cbcb + c + b));

            const scalar second_part =
                T * dVdT_p(p, T, z) * dd2p_dT2_V(p, T, z);
            const scalar dCvR_dT = a_log_part + second_part;

            const scalar d_dT_p__dpdT_V =
                dd2p_dT2_V(p, T, z) + dd2p_dVdT(p, T, z) * dVdT_p(p, T, z);

            const scalar dCpRMCvR_dT =
                dVdT_p(p, T, z) * (dpdT_V(p, T, z) + T * d_dT_p__dpdT_V) +
                T * dpdT_V(p, T, z) * dd2V_dT2_p(p, T, z);

            const scalar dCpR_dT = dCvR_dT + dCpRMCvR_dT;

            return Foam::residualSpecificHeatDerivative{
                .dCpdT_R = dCpR_dT * W,
                .dCvdT_R = dCvR_dT * W,
                .dCp_minus_Cv_dT_R = dCpRMCvR_dT * W,
                .Z = z,
            };
        }

        [[nodiscard]] inline scalar e(scalar p, scalar T, scalar W) const
        {
            return this->residualHE(p, T, W).e_R;
        }

        [[nodiscard]] inline scalar h(scalar p, scalar T, scalar W) const
        {
            return this->residualHE(p, T, W).h_R;
        }

        //- Return Cp contribution [J/(kg K)]
        inline scalar Cp(scalar p, scalar T, scalar W) const
        {
            return this->residualSpecificHeat(p, T, W).Cp_R;
        }

        //- Return Cv contribution [J/(kg K)]
        inline scalar Cv(scalar p, scalar T, scalar W) const
        {
            return this->residualSpecificHeat(p, T, W).Cv_R;
        }

        //- Return dCv/dT_p contribution [J/(kg K^2)]
        inline scalar dCvdT(scalar p, scalar T, scalar W) const
        {
            return this->residualSpecificHeatDerivative(p, T, W).dCvdT_R;
        }

        //- Return dCp/dT_p contribution [J/(kg K^2)]
        inline scalar dCpdT(scalar p, scalar T, scalar W) const
        {
            return this->residualSpecificHeatDerivative(p, T, W).dCpdT_R;
        }

        //- Return entropy contribution to the integral of Cp/T [J/kg/K]
        inline scalar sp(const scalar p, const scalar T, scalar W) const
        {
            return 0;
        }

        //- Return entropy contribution to the integral of Cv/T [J/kg/K]
        inline scalar sv(const scalar p, const scalar T, scalar W) const
        {
            return 0;
        }

        //- Return compressibility [s^2/m^2]
        inline scalar psi(scalar p, scalar T, scalar W) const
        {
            return 1 / (this->Z(p, T) * this->R(W) * T);
        }

        //- Return (Cp - Cv) [J/(kg K]
        inline scalar CpMCv(scalar p, scalar T, scalar W) const
        {
            // CpMCv should return Cp-Cv of real gas, instead of CpR-CvR.
            return this->residualSpecificHeat(p, T, W).Cp_R_minus_Cv_R +
                this->R(W);
        }

        inline scalar alphav(scalar p, scalar T, scalar W) const
        {
            const scalar Z = this->Z(this->dimlessCoeff(p, T));
            const scalar V = Z * RR * T / p;
            const scalar dVdT_const_p = this->dVdT_p(p, T, Z);

            const scalar alpha_v = dVdT_const_p / V;
            return alpha_v;
        }
        //------------------ Derivatives -------------------------------------//

        [[nodiscard]] inline scalar dpdT_V(scalar p, scalar T, scalar z) const
        {
            const scalar V = z * RR * T / p;
            const scalar result =
                RR / (V - b) - da_dT / (V * (V + b) + c * (V - b));
            return result;
        }
        [[nodiscard]] inline scalar dTdp_V(scalar p, scalar T, scalar Z) const
        {
            return 1 / dpdT_V(p, T, Z);
        }

        [[nodiscard]] inline scalar dpdV_T(scalar p, scalar T, scalar z) const
        {
            const scalar V = z * RR * T / p;
            const scalar result = -RR * T / sqr(V - b) +
                a * (2 * V + b + c) / sqr(sqr(V) + (b + c) * V - b * c);
            return result;
        }
        [[nodiscard]] inline scalar dVdp_T(scalar p, scalar T, scalar z) const
        {
            return 1 / dpdV_T(p, T, z);
        }

        [[nodiscard]] inline scalar dVdT_p(scalar p, scalar T, scalar z) const
        {
            return -1 * dpdT_V(p, T, z) * dVdp_T(p, T, z);
        }
        [[nodiscard]] inline scalar dTdV_p(scalar p, scalar T, scalar z) const
        {
            return 1 / dpdT_V(p, T, z);
        }

        [[nodiscard]] inline scalar dd2p_dT2_V(scalar p, scalar T,
                                               scalar z) const
        {
            const scalar V = z * RR * T / p;
            const scalar result = -dd2a_dT2 / (V * (V + b) + c * (V - b));
            return result;
        }

        [[nodiscard]] inline scalar dd2V_dT2_p(scalar p, scalar T,
                                               scalar z) const
        {
            const scalar V = z * RR * T / p;
            const scalar Vp = dVdT_p(p, T, z);
            // denominator
            const scalar D = 3 * p * sqr(V) + (2 * c * p - 2 * RR * T) * V + a +
                (-RR * c - RR * b) * T + (-2 * b * c - sqr(b)) * p;

            const scalar numerator1 = 2 * RR * V * Vp +
                (-da_dT + RR * (c + b)) * Vp - dd2a_dT2 * V + b * dd2a_dT2;

            const scalar term1 = numerator1 / D;

            const scalar factor1 = RR * sqr(V) + (-da_dT + RR * (c + b)) * V +
                b * da_dT - RR * b * c;

            const scalar factor2 = 6 * p * V * Vp +
                (2 * c * p - 2 * RR * T) * Vp - 2 * RR * V + da_dT - RR * c -
                RR * b;

            const scalar term2 = (factor1 * factor2) / sqr(D);

            const scalar result = term1 - term2;
            return result;
        }

        [[nodiscard]] inline scalar dd2p_dVdT(scalar p, scalar T,
                                              scalar z) const
        {
            const scalar V = z * RR * T / p;
            const scalar left_part = -RR / sqr(V - b);
            const scalar a_part =
                da_dT * (2 * V + b + c) / sqr(sqr(V) + (b + c) * V - b * c);
            const scalar result = left_part + a_part;
            return result;
        }

        [[nodiscard]] inline scalar dhdp_T(scalar p, scalar T, scalar W) const
        {
            const scalar z = this->Z(p, T);
            const scalar V = z * RR * T / p;
            const scalar dHdp_T = V - T * dVdT_p(p, T, z);
            return dHdp_T / this->R(W);
        }

        [[nodiscard]] inline scalar dedp_T(scalar p, scalar T, scalar W) const
        {
            const scalar z = this->Z(p, T);
            [[maybe_unused]] const scalar V = z * RR * T / p;
            const scalar dUdp_T = -T * dVdT_p(p, T, z) - p * dVdp_T(p, T, z);
            return dUdp_T / this->R(W);
        }
    };

  public:
    inline PatelTejaValderrama(const Specie &sp, const realFluidProperty &prop)
        : Specie{sp}, property_{prop}
    {}
    PatelTejaValderrama(const PatelTejaValderrama &) = default;
    // May be not required and not implemented?
    //    PatelTejaValderrama(const Specie &sp);

    inline PatelTejaValderrama(const word &name, const dictionary &dict)
        : Specie{name, dict}, property_{dict}
    {
        property_.requireRealGasEOS(name, true);
    }

    inline PatelTejaValderrama(const word &name, const PatelTejaValderrama &ptv)
        : Specie{name, ptv}, property_{ptv.property_}
    {}

    //- Construct and return a clone
    inline autoPtr<PatelTejaValderrama> clone() const
    {
        return autoPtr<PatelTejaValderrama<Specie>>(
            new PatelTejaValderrama<Specie>(*this));
    }

    inline scalar Tc() const noexcept { return this->property_.Tc_; }
    inline scalar Pc() const noexcept { return this->property_.Pc_; }
    inline scalar Vc() const noexcept { return this->property_.Vc_; }
    inline scalar Zc() const noexcept { return this->property_.Zc(); }
    inline scalar omega() const noexcept { return this->property_.omega_; }
    inline bool preferGas() const noexcept { return this->property_.is_gas; }

    inline static word typeName()
    {
        return "PatelTejaValderrama<" + word(Specie::typeName_()) + '>';
    }

    inline EOSCore core(scalar p, scalar T) const
    {
        const scalar Tr = T / Tc();
        assert(Tr > 0);
        const scalar Omega_a = 0.66121 - 0.761057 * Zc();
        const scalar Omega_b = 0.02207 + 0.20868 * Zc();
        const scalar Omega_c = 0.57765 - 1.87080 * Zc();
        const scalar F = 0.46283 + 3.58230 * omega() * Zc() +
            8.19417 * sqr(omega() * sqr(Zc()));

        const scalar a_c = Omega_a * sqr(RR * Tc()) / Pc();
        const scalar alpha_Tr = sqr(1 + F * (1 - sqrt(Tr)));
        const scalar a = a_c * alpha_Tr;
        const scalar dalpha_dTr = F * (F * (sqrt(Tr) - 1) - 1) / sqrt(Tr);
        const scalar dd2alpha_dTr2 = F * (F + 1) / (2 * Tr * sqrt(Tr));
        const scalar da_dT = a_c / Tc() * dalpha_dTr;
        const scalar dd2a_dT2 = a_c / sqr(Tc()) * dd2alpha_dTr2;

        const scalar b = Omega_b * RR * Tc() / Pc();
        const scalar c = Omega_c * RR * Tc() / Pc();
        assert(b > 0);

        return EOSCore{
            .a = a,
            .b = b,
            .c = c,
            .da_dT = da_dT,
            .dd2a_dT2 = dd2a_dT2,
            .is_gas = property_.is_gas,
        };
    }

    //- Is the equation of state is incompressible i.e. rho != f(p)
    static constexpr bool incompressible = false;

    //- Is the equation of state is isochoric i.e. rho = const
    static constexpr bool isochoric = false;

    //- Return compression factor []
    inline scalar Z(scalar p, scalar T) const
    {
        return this->core(p, T).Z(p, T);
    }

    //- Return density [kg/m^3]
    inline scalar rho(scalar p, scalar T) const
    {
        return this->core(p, T).rho(p, T, this->W());
    }

    //- Return enthalpy contribution [J/kg]
    inline scalar h(const scalar p, const scalar T) const
    {
        return this->core(p, T).h(p, T, this->W());
    }

    //- Return Cp contribution [J/(kg K)]
    inline scalar Cp(scalar p, scalar T) const
    {
        return this->core(p, T).Cp(p, T, this->W());
    }

    //- Return dCp/dT_p contribution [J/(kg K^2)]
    inline scalar dCpdT(const scalar p, const scalar T) const
    {
        return this->core(p, T).dCpdT(p, T, this->W());
    }

    //- Return internal energy contribution [J/kg]
    inline scalar e(const scalar p, const scalar T) const
    {
        return this->core(p, T).e(p, T, this->W());
    }

    //- Return dCv/dT_p contribution [J/(kg K^2)]
    inline scalar dCvdT(const scalar p, const scalar T) const
    {
        return this->core(p, T).dCvdT(p, T, this->W());
    }

    //- Return Cv contribution [J/(kg K)]
    inline scalar Cv(scalar p, scalar T) const
    {
        return this->core(p, T).Cv(p, T, this->W());
    }

    //- Return entropy contribution to the integral of Cp/T [J/kg/K]
    inline scalar sp(const scalar p, const scalar T) const
    {
        return this->core(p, T).sp(p, T, this->W());
    }

    //- Return entropy contribution to the integral of Cv/T [J/kg/K]
    inline scalar sv(const scalar p, const scalar T) const
    {
        return this->core(p, T).sv(p, T, this->W());
    }

    //- Return compressibility [s^2/m^2]
    inline scalar psi(scalar p, scalar T) const
    {
        return this->core(p, T).psi(p, T, this->W());
    }

    //- Return (Cp - Cv) [J/(kg K]
    inline scalar CpMCv(scalar p, scalar T) const
    {
        return this->core(p, T).CpMCv(p, T, this->W());
    }

    //- Return volumetric coefficient of thermal expansion [1/T]
    inline scalar alphav(const scalar p, const scalar T) const
    {
        return this->core(p, T).alphav(p, T, this->W());
    }

    inline scalar dhdp_T(const scalar p, const scalar T) const
    {
        return this->core(p, T).dhdp_T(p, T, this->W());
    }

    inline scalar dedp_T(const scalar p, const scalar T) const
    {
        return this->core(p, T).dedp_T(p, T, this->W());
    }

    //- Write to Ostream
    void write(Ostream &os) const;

    // Member Operators

    void operator+=(const PatelTejaValderrama &);

    void operator*=(const scalar);

    // Friend operators

    friend PatelTejaValderrama operator+
        <Specie>(const PatelTejaValderrama &, const PatelTejaValderrama &);

    friend PatelTejaValderrama operator*
        <Specie>(const scalar s, const PatelTejaValderrama &);

    friend PatelTejaValderrama operator==
        <Specie>(const PatelTejaValderrama &, const PatelTejaValderrama &);

    // Ostream Operator

    friend Ostream &operator<< <Specie>(Ostream &, const PatelTejaValderrama &);
};

} // namespace Foam

template <class Specie>
inline void Foam::PatelTejaValderrama<Specie>::operator+=(
    const PatelTejaValderrama &)
{
    noCoefficientMixing(PatelTejaValderrama);
}

template <class Specie>
inline void Foam::PatelTejaValderrama<Specie>::operator*=(const scalar)
{
    noCoefficientMixing(PatelTejaValderrama);
}

template <class Specie>
inline Foam::PatelTejaValderrama<Specie> Foam::operator+(
    const Foam::PatelTejaValderrama<Specie> &,
    const Foam::PatelTejaValderrama<Specie> &)
{
    noCoefficientMixing(PatelTejaValderrama);
}

template <class Specie>
inline Foam::PatelTejaValderrama<Specie> Foam::operator*(
    const Foam::scalar s, const Foam::PatelTejaValderrama<Specie> &)
{
    noCoefficientMixing(PatelTejaValderrama);
}

template <class Specie>
inline Foam::PatelTejaValderrama<Specie> Foam::operator==(
    const Foam::PatelTejaValderrama<Specie> &,
    const Foam::PatelTejaValderrama<Specie> &)
{
    noCoefficientMixing(PatelTejaValderrama);
}

#ifdef NoRepository
#include "PatelTejaValderrama.C"
#endif

#endif // MIKENO_PATELTEJAVALDERRAMA_H
