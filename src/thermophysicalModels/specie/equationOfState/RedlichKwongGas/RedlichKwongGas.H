/*
// Created by joseph on 25-7-11.
//

License
    This file is part of Mikeno.

Class
    Foam::RedlichKwongGas

Description
    RedlichKwongGas cubic equation of state for gases

    Coefficient mixing is very inaccurate and not supported, so this equation of
    state can only be applied to mixtures by using a thermophysical properties
    model which mixes specie property value.

Usage
    \table
        Property     | Description
        Tc           | Critical temperature [K]
        Vc           | Critical volume [m^3/kmol]
        Pc           | Critical pressure [Pa]
        phase        | Preferred phase (vapor/liquid)
    \endtable

    Example specification of the RedlichKwongGas equation of state for methane:
    \verbatim
    equationOfState
    {
        Tc        190.55;
        Vc        0.0285;
        Pc        4.595e6;
        phase     vapor;
    }
    \endverbatim

SourceFiles
    RedlichKwongGasI.H
    RedlichKwongGas.C

\*---------------------------------------------------------------------------*/

#ifndef RedlichKwongGas_H
#define RedlichKwongGas_H

#include "autoPtr.H"
#include "GeometricSymmTensorField.H"
#include "../equationOfState/cubicEOS.H"

namespace Foam {

template<class Specie> class RedlichKwongGas;

template<class Specie>
inline RedlichKwongGas<Specie> operator+
(
        const RedlichKwongGas<Specie>&,
        const RedlichKwongGas<Specie>&
);

template<class Specie>
inline RedlichKwongGas<Specie> operator*
(
        const scalar,
        const RedlichKwongGas<Specie>&
);

template<class Specie>
inline RedlichKwongGas<Specie> operator==
(
        const RedlichKwongGas<Specie>&,
        const RedlichKwongGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
        Ostream&,
        const RedlichKwongGas<Specie>&
);

template<class Specie>
class RedlichKwongGas : public Specie {

public:

    realFluidProperty property_;

    struct EOSCore {
        scalar a;
        scalar b;
        bool is_gas;

    protected:
        //- Molar volume [m^3/kmol]
    public:
        inline cubicEOSDimlessCoefficient dimlessCoeff(scalar p, scalar T) const {
            const scalar A = a *p/ sqr(RR*T)/ sqrt(T);
            const scalar B = b * p / (RR * T);
            return {A, B};
        }
        // compute from AB, that is real EOSCore
        inline scalar Z(const cubicEOSDimlessCoefficient&AB) const {
            const scalar A = AB.A;
            const scalar B = AB.B;

            const scalar a2 = -1;
            const scalar a1 = A - B - sqr(B);
            const scalar a0 = -A * B;

            const scalar Z= solveCubicEquation(a2, a1, a0, this->is_gas, B);
            assert(Z>B);
            return Z;
        }
        //- Return residual h and e
        inline Foam::residualEnergy residualHE(const scalar p, const scalar T, scalar W) const {
            const auto coeff = this->dimlessCoeff(p, T);

            const scalar A = coeff.A;
            const scalar B = coeff.B;
            const scalar Z = this->Z(coeff);

            const scalar e_R_dimless = -(3 * A) / (2 * B) * log(Z / (Z + B));
            const scalar h_R_dimless =
                    (Z - 1)
                    + e_R_dimless;
            const scalar e_R = T * this->R(W) * e_R_dimless;
            const scalar h_R = T * this->R(W) * h_R_dimless;
            return Foam::residualEnergy{h_R, e_R, Z};
        }

        //- Return residual Cp and Cv
        inline Foam::residualSpecificHeat residualSpecificHeat(const scalar p, const scalar T, scalar W) const {
            const auto coeff = this->dimlessCoeff(p, T);

            const scalar A = coeff.A;
            const scalar B = coeff.B;
            const scalar Z = this->Z(coeff);
            const scalar Cv_dimless = (3 * A) / (4 * B) * log(1 + B / Z);

            const scalar CpRMCvR_dimless =
                    -1 -
                    sqr(1 / (Z - B) + A / (2 * Z * (Z + B)))
                    /
                    (-1 / sqr(Z - B) + A * (2 * Z + B) / sqr(Z) / sqr(Z + B));
            const scalar Cp_dimless = Cv_dimless + CpRMCvR_dimless;

            return Foam::residualSpecificHeat{
                    Cp_dimless * this->R(W),
                    Cv_dimless * this->R(W),
                    CpRMCvR_dimless * this->R(W),
                    Z
            };
        }

        inline Foam::residualSpecificHeatDerivative residualSpecificHeatDerivative(const scalar p, const scalar T, const scalar W) const {
            const auto coeff = this->dimlessCoeff(p, T);
            const scalar z = this->Z(coeff);
            const scalar V = z*RR*T/p;

            const scalar dCvdT_p = (3*a)/(4*b*sqr(T)*sqrt(T)) * log(1+b/V) + T* dVdT_p(p,T,z) * dd2p_dT2_V(p,T,z);
            #warning "TODO: implement derivative of residual specific heat difference d(Cp^R-Cv^R)/dT. Ignored currently"
            const scalar dCpdT_p = dCvdT_p + 0;

            return Foam::residualSpecificHeatDerivative{
                .dCpdT_R = dCpdT_p * W,
                .dCvdT_R = dCvdT_p * W,
                .dCp_minus_Cv_dT_R = 0 * W,
                .Z = z,
            };
        }

        inline scalar R(scalar W) const {
            return RR/W;
        }
        // compute from p and T
        inline scalar Z(scalar p, scalar T) const {
            auto AB=this->dimlessCoeff(p,T);
            return this->Z(AB);
        }

        inline scalar rho(scalar p, scalar T, scalar W) const {
            const scalar Z=this->Z(p,T);
            return p / (Z * this->R(W) * T);
        }

        //- Return enthalpy contribution [J/kg]
        inline scalar h(const scalar p, const scalar T, scalar W) const {
            return this->residualHE(p,T,W).h_R;
        }
        //- Return internal energy contribution [J/kg]
        inline scalar e(const scalar p, const scalar T, scalar W) const {
            return this->residualHE(p,T,W).e_R;
        }

        //- Return Cp contribution [J/(kg K]
        inline scalar Cp(scalar p, scalar T, scalar W) const {
            return this->residualSpecificHeat(p,T,W).Cp_R;
        }


        //- Return Cv contribution [J/(kg K]
        inline scalar Cv(scalar p, scalar T, scalar W) const{
            return this->residualSpecificHeat(p,T,W).Cv_R;
        }

        //- Return entropy contribution to the integral of Cp/T [J/kg/K]
        inline scalar sp(const scalar p, const scalar T, scalar W) const {
            return 0;
        }

        //- Return entropy contribution to the integral of Cv/T [J/kg/K]
        inline scalar sv(const scalar p, const scalar T, scalar W) const {
            return 0;
        }

        //- Return compressibility [s^2/m^2]
        inline scalar psi(scalar p, scalar T, scalar W) const {
            return 1/(this->Z(p,T)*this->R(W)*T);
        }

        //- Return (Cp - Cv) [J/(kg K]
        inline scalar CpMCv(scalar p, scalar T, scalar W) const{
            // CpMCv should return Cp-Cv of real gas, instead of CpR-CvR.
            return this->residualSpecificHeat(p,T,W).Cp_R_minus_Cv_R + this->R(W);
        }

        //- Return dCv/dT_p contribution [J/(kg K^2)]
        inline scalar dCvdT(scalar p, scalar T, scalar W) const {
            return this->residualSpecificHeatDerivative(p,T,W).dCvdT_R;
        }

        //- Return dCp/dT_p contribution [J/(kg K^2)]
        inline scalar dCpdT(scalar p, scalar T, scalar W) const {
            return this->residualSpecificHeatDerivative(p,T,W).dCpdT_R;
        }

        //- Return volumetric coefficient of thermal expansion [1/T]
        inline scalar alphav(const scalar p, const scalar T) const;

        //- Return isobaric expansion coeff [m^3/mol/K]
        inline scalar dVdT_p(const scalar p, const scalar T, const scalar z) const {
            return - dpdT_V(p, T, z) / dpdV_T(p, T, z);
        }

        inline scalar dTdV_p(const scalar p, const scalar T, const scalar z) const {
            return 1./dVdT_p(p, T, z);
        }

        //- Return isothermal compress coeff [m^3/mol/Pa]
        inline scalar dVdp_T(const scalar p, const scalar T, const scalar z) const {
            return 1./dpdV_T(p, T, z);
        }

        inline scalar dpdV_T(const scalar p, const scalar T, const scalar z) const {
            // return 1./dVdp_T(p, T, z);
            const scalar V = z*RR*T/p;
            const scalar result = -RR*T/sqr(V-b) + a*(2*V+b)/(sqrt(T)*sqr(V*(V+b)));
            return result;
        }

        inline scalar dpdT_V(const scalar p, const scalar T, const scalar z) const {
            const scalar V = z*RR*T/p;
            const scalar result = RR/(V-b) + a/(2*T*sqrt(T) * V*(V+b));
            return result;
        }

        inline scalar dTdp_V(const scalar p, const scalar T, const scalar z) const {
            return 1./dpdT_V(p, T, z);
        }

        inline scalar dd2p_dT2_V(const scalar p, const scalar T, const scalar z) const {
            const scalar V = z*RR*T/p;
            const scalar result = - 3*a / (4*sqr(T) * sqrt(T) * V*(V+b));
            return result;
        }

        //- Return isothermal enthalpy-pressure derivative [m/kg]
        inline scalar dhdp_T(const scalar p, const scalar T, scalar W) const {
            const auto z=this->Z(p, T);
            const scalar V = z*RR*T/p;
            const scalar dHdp_T = V - T*dVdT_p(p, T, z);
            return dHdp_T/W;
        }

        //- Return isothermal internal energy-pressure derivative [m/kg]
        inline scalar dedp_T(const scalar p, const scalar T, scalar W) const {
            const auto z=this->Z(p, T);
            [[maybe_unused]] const scalar V = z*RR*T/p;
            const scalar dUdp_T = - T*dVdT_p(p, T, z) - p*dVdp_T(p, T, z);
            return dUdp_T/W;
        }
    };

public:
    inline RedlichKwongGas(
            const Specie &sp,
            const realFluidProperty&prop
       );

    RedlichKwongGas(const word&name, const dictionary&dict);

    inline RedlichKwongGas(const word&name, const RedlichKwongGas&);

    inline autoPtr<RedlichKwongGas> clone() const;


    // Member Functions
    inline scalar Tc() const noexcept {
        return this->property_.Tc_;
    }
    inline scalar Pc() const noexcept {
        return this->property_.Pc_;
    }
    inline scalar Vc() const noexcept {
        return this->property_.Vc_;
    }
    inline scalar Zc() const noexcept {
        return this->property_.Zc();
    }
    inline bool preferGas() const noexcept {
        return this->property_.is_gas;
    }

    //- Return the instantiated type name
    static word typeName()
    {
        return "RedlichKwongGas<" + word(Specie::typeName_()) + '>';
    }
//    inline Foam::cubicEOSRawCoefficient coeffRaw(scalar p, scalar T) const;

    inline EOSCore core(scalar p[[maybe_unused]], scalar T[[maybe_unused]]) const {
        const auto Tc_ = this->Tc();
//    const auto Vc_=this->Vc_;
        const auto Pc_ = this->Pc();
        const scalar a = 0.42748025 * sqr(RR * Tc_) * sqrt(Tc_) / Pc_;
        const scalar b = 0.08664035 * RR * Tc_ / Pc_;
        return {a, b, this->preferGas()};
    }

    //- Is the equation of state is incompressible i.e. rho != f(p)
    static const bool incompressible = false;

    //- Is the equation of state is isochoric i.e. rho = const
    static const bool isochoric = false;

    //- Return compression factor []
    inline scalar Z(scalar p, scalar T) const{
        return this->core(p,T).Z(p,T);
    }

    //- Return density [kg/m^3]
    inline scalar rho(scalar p, scalar T) const {
        return this->core(p,T).rho(p,T,this->W());
    }

    //- Return enthalpy contribution [J/kg]
    inline scalar h(const scalar p, const scalar T) const {
        return this->core(p,T).h(p,T,this->W());
    }

    //- Return Cp contribution [J/(kg K)]
    inline scalar Cp(scalar p, scalar T) const{
        return this->core(p,T).Cp(p,T,this->W());
    }

    //- Return dCp/dT_p contribution [J/(kg K^2)]
    inline scalar dCpdT(const scalar p, const scalar T) const {
        return this->core(p,T).dCpdT(p,T,this->W());
    }

    //- Return internal energy contribution [J/kg]
    inline scalar e(const scalar p, const scalar T) const{
        return this->core(p,T).e(p,T,this->W());
    }

    //- Return Cv contribution [J/(kg K)]
    inline scalar Cv(scalar p, scalar T) const{
        return this->core(p,T).Cv(p,T,this->W());
    }

    //- Return dCv/dT_p contribution [J/(kg K^2)]
    inline scalar dCvdT(const scalar p, const scalar T) const {
        return this->core(p,T).dCvdT(p,T,this->W());
    }

    //- Return entropy contribution to the integral of Cp/T [J/kg/K]
    inline scalar sp(const scalar p, const scalar T) const {
        return this->core(p,T).sp(p,T,this->W());
    }

    //- Return entropy contribution to the integral of Cv/T [J/kg/K]
    inline scalar sv(const scalar p, const scalar T) const{
        return this->core(p,T).sv(p,T,this->W());
    }

    //- Return compressibility [s^2/m^2]
    inline scalar psi(scalar p, scalar T) const{
        return this->core(p,T).psi(p,T,this->W());
    }

    //- Return (Cp - Cv) [J/(kg K]
    inline scalar CpMCv(scalar p, scalar T) const{
        return this->core(p,T).CpMCv(p,T,this->W());
    }

    //- Return volumetric coefficient of thermal expansion [1/T]
    inline scalar alphav(const scalar p, const scalar T) const{
        return this->core(p,T).alphav(p,T,this->W());
    }

    //- Return isothermal enthalpy-pressure derivative [m/kg]
    inline scalar dhdp_T(const scalar p, const scalar T) const {
        return this->core(p, T).dhdp_T(p, T, this->W());
    }

    //- Return isothermal internal energy-pressure derivative [m/kg]
    inline scalar dedp_T(const scalar p, const scalar T) const {
        return this->core(p, T).dedp_T(p, T, this->W());
    }


    class EOSMixer {
    public:
        EOSMixer(const speciesTable &, const dictionary&) { }

        template<class getCoreFun, class getXiFun>
        inline EOSCore mix(label nSpecie,
                           getCoreFun &coreFun,
                           getXiFun &specieMoleFracFun) const {
            scalar sqrt_a_sum = 0;
            scalar bSum = 0;
            const bool is_gas=coreFun(0).is_gas;
            for (label spIdx = 0; spIdx < nSpecie; spIdx++) {
                const EOSCore core_i = coreFun(spIdx);
                const scalar Xi = specieMoleFracFun(spIdx);
                sqrt_a_sum += Xi * sqrt(core_i.a);
                bSum += Xi * core_i.b;
            }
            return EOSCore{sqr(sqrt_a_sum), bSum, is_gas};
        }
    };

    // IO

    //- Write to Ostream
    void write(Ostream& os) const;

    // Member Operators

    inline void operator+=(const RedlichKwongGas&);
    inline void operator*=(const scalar);


    // Friend operators

    friend RedlichKwongGas operator+ <Specie>
            (
                    const RedlichKwongGas&,
                    const RedlichKwongGas&
            );

    friend RedlichKwongGas operator* <Specie>
            (
                    const scalar s,
                    const RedlichKwongGas&
            );

    friend RedlichKwongGas operator== <Specie>
            (
                    const RedlichKwongGas&,
                    const RedlichKwongGas&
            );


    // Ostream Operator

    friend Ostream& operator<< <Specie>
            (
                    Ostream&,
                    const RedlichKwongGas&
            );
};

} // End namespace Foam


#include "RedlichKwongGasI.H"

#ifdef NoRepository
#include "RedlichKwongGas.C"
#endif

#endif //RedlichKwongGas_H
