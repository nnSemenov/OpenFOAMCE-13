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

Class
    Foam::PengRobinsonGas

Description
    PengRobinsonGas cubic equation of state for gases

    Coefficient mixing is very inaccurate and not supported, so this equation of
    state can only be applied to mixtures by using a thermophysical properties
    model which mixes specie property value.

Usage
    \table
        Property     | Description
        Tc           | Critical temperature [K]
        Vc           | Critical volume [m^3/kmol]
        Pc           | Critical pressure [Pa]
        omega        | Acentric factor [-]
        phase        | Preferred phase (vapor/liquid)
    \endtable

    Example specification of the PengRobinsonGas equation of state for methane:
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
    PengRobinsonGasI.H
    PengRobinsonGas.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonGas_H
#define PengRobinsonGas_H

#include "autoPtr.H"
#include "speciesTable.H"
#include "../equationOfState/cubicEOS.H"
#include <cassert>
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// Forward declaration of friend functions and operators

    template<class Specie>
    class PengRobinsonGas;

    template<class Specie>
    inline PengRobinsonGas<Specie> operator+
            (
                    const PengRobinsonGas<Specie> &,
                    const PengRobinsonGas<Specie> &
            );

    template<class Specie>
    inline PengRobinsonGas<Specie> operator*
            (
                    const scalar,
                    const PengRobinsonGas<Specie> &
            );

    template<class Specie>
    inline PengRobinsonGas<Specie> operator==
            (
                    const PengRobinsonGas<Specie> &,
                    const PengRobinsonGas<Specie> &
            );

    template<class Specie>
    Ostream &operator<<
            (
                    Ostream &,
                    const PengRobinsonGas<Specie> &
            );


/*---------------------------------------------------------------------------*\
                       Class PengRobinsonGas Declaration
\*---------------------------------------------------------------------------*/

    template<class Specie>
    class PengRobinsonGas
            :
                    public Specie {
        // Private Data
        realFluidProperty property_;

    public:

        struct EOSCore {
            scalar a;
            scalar b;
            scalar da_dT;
            scalar dd2a_dT2;
            bool is_gas;

            inline cubicEOSDimlessCoefficient dimlessCoeff(scalar p, scalar T) const {
                const scalar A = a * p / sqr(RR * T);
                const scalar B = b * p / (RR * T);
                return {A,B};
            }

            inline scalar R(scalar W) const {
                return RR/W;
            }

            inline scalar Z(const cubicEOSDimlessCoefficient&AB) const {
                const scalar A = AB.A;
                const scalar B = AB.B;

                const scalar a2 = B-1;
                const scalar a1=A -2*B -3*sqr(B);
                const scalar a0 = pow3(B) + sqr(B)-A*B;
                const scalar Z = solveCubicEquation(a2,a1,a0,this->is_gas, B);
                assert(Z>B);
                return Z;
            }

            inline scalar Z(scalar p, scalar T) const {
                auto AB=this->dimlessCoeff(p,T);
                return this->Z(AB);
            }

            inline scalar rho(scalar p, scalar T, scalar W) const {
                const scalar Z=this->Z(p,T);
                return p / (Z * this->R(W) * T);
            }

            inline Foam::residualEnergy residualHE(scalar p, scalar T, scalar W) const {
                const auto AB=this->dimlessCoeff(p,T);
                [[maybe_unused]] const scalar A=AB.A;
                const scalar B=AB.B;
                const scalar Z=this->Z(AB);

                const scalar log_part = log(
                            (Z-(sqrt(2.0)-1)*B)
                            /(Z+(sqrt(2.0)+1)*B)
                        );
                const scalar aT_part = a - T*this->da_dT;
                const scalar e_R_dimless = 1/(2* sqrt(2.0)*b*RR*T)*aT_part*log_part;
                const scalar h_R_dimless = Z-1 + e_R_dimless;
                return Foam::residualEnergy{
                        T * this->R(W) * h_R_dimless,
                        T * this->R(W) * e_R_dimless,
                        Z
                };
            }

            inline Foam::residualSpecificHeat residualSpecificHeat(scalar p, scalar T, scalar W) const {
                const auto AB=this->dimlessCoeff(p,T);
                const scalar B=AB.B;
                const scalar Z=this->Z(AB);

                const scalar log_part = log(
                        (Z-(sqrt(2.0)-1)*B)
                        /(Z+(sqrt(2.0)+1)*B)
                );
                const scalar V = Z*RR*T/p;
                [[maybe_unused]] const scalar VVBbbVb=V*(V+b)+b*(V-b);
                const scalar pB_part = p/(2*sqrt(2.0)*sqr(RR)*B);

                const scalar Cv_dimless = -this->dd2a_dT2*pB_part*log_part;

                // compute CpMCv
                const scalar CpMCv_dimless = -1 +T/RR * dpdT_V(p, T, Z) * dVdT_p(p, T, Z);
                const scalar Cp_dimless = Cv_dimless+CpMCv_dimless;

                return Foam::residualSpecificHeat{
                        Cp_dimless * this->R(W),
                        Cv_dimless * this->R(W),
                        CpMCv_dimless * this->R(W),
                        Z
                };
            }

            inline Foam::residualSpecificHeatDerivative residualSpecificHeatDerivative(scalar p, scalar T, scalar W) const {
                const scalar Z=this->Z(this->dimlessCoeff(p,T));

                const scalar V = Z*RR*T/p;
                const scalar log_part = std::log(
                (V-(sqrt(2.0)-1)*b)
                /(V+(sqrt(2.0)+1)*b)
                );
                // J/kmol/K
                const scalar dCvR_dT_first_part = -dd2a_dT2 /(2*sqrt(2.0)*b) * log_part;
                const scalar dCvR_dT_second_part = T*dVdT_p(p, T, Z)*dd2p_dT2_V(p, T, Z);
                const scalar dCvR_dT = dCvR_dT_first_part + dCvR_dT_second_part;
                #warning "TODO: implement derivative of residual specific heat difference d(Cp^R-Cv^R)/dT. Ignored currently"
                const scalar dCpR_dT = dCvR_dT + 0;
                return Foam::residualSpecificHeatDerivative{
                    .dCpdT_R = dCpR_dT * W,
                    .dCvdT_R = dCvR_dT * W,
                    .dCp_minus_Cv_dT_R = 0 * W,
                    .Z=Z,
                };
            }

            inline scalar e(scalar p, scalar T, scalar W) const {
                return this->residualHE(p,T,W).e_R;
            }

            inline scalar h(scalar p, scalar T, scalar W) const {
                return this->residualHE(p,T,W).h_R;
            }

            //- Return Cp contribution [J/(kg K)]
            inline scalar Cp(scalar p, scalar T, scalar W) const {
                return this->residualSpecificHeat(p,T,W).Cp_R;
            }

            //- Return Cv contribution [J/(kg K)]
            inline scalar Cv(scalar p, scalar T, scalar W) const{
                return this->residualSpecificHeat(p,T,W).Cv_R;
            }

            //- Return dCv/dT_p contribution [J/(kg K^2)]
            inline scalar dCvdT(scalar p, scalar T, scalar W) const {
                return this->residualSpecificHeatDerivative(p,T,W).dCvdT_R;
            }

            //- Return dCp/dT_p contribution [J/(kg K^2)]
            inline scalar dCpdT(scalar p, scalar T, scalar W) const {
                return this->residualSpecificHeatDerivative(p,T,W).dCpdT_R;
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

            inline scalar alphav(const scalar p, const scalar T) const;

            inline scalar dpdT_V(const scalar p, const scalar T,
                                 const scalar z) const {
                const scalar V = z*RR*T/p;
                const scalar result = RR/(V-b) - da_dT/(V*(V+b)+b*(V-b));
                return result;
            }

            inline scalar dTdp_V(const scalar p, const scalar T, const scalar Z) const {
                return 1./dpdT_V(p,T,Z);
            }

            inline scalar dVdT_p(const scalar p, const scalar T,
                                 const scalar z) const {
              const scalar V = z * RR * T / p;
              const scalar dVdT_p = (RR * sqr(V) + (2 * RR * b - da_dT) * V +
                                     b * (da_dT - RR * b)) /
                                    (3 * p * sqr(V) + 2 * V * (p * b - RR * T) +
                                     a - 2 * b * RR * T - 3 * p * sqr(b));
              return dVdT_p;
            }

            inline scalar dTdV_p(const scalar p, const scalar T, const scalar Z) const {
                return 1./dVdT_p(p,T,Z);
            }

            inline scalar dVdp_T(const scalar p, const scalar T,
                                 const scalar z) const {
              const scalar V = z * RR * T / p;
              const scalar dVdp_T =
                  (-RR * T * sqr(V) - (2 * RR * T * b - a) * V - a * b +
                   sqr(b) * RR * T) /
                  (3 * sqr(p * V) + 2 * p * V * (p * b - RR * T) + a * p -
                   2 * b * RR * T * p - 3 * sqr(p * b));
              return dVdp_T;
            }

            inline scalar dpdV_T(const scalar p, const scalar T,
                                 const scalar z) const {
                return 1./dVdp_T(p, T, z);
            }

            inline scalar dd2p_dT2_V(const scalar p, const scalar T,
                                 const scalar z) const {
                const scalar V = z * RR * T / p;
                const scalar result = -this->dd2a_dT2 / (V*(V+b)+b*(V-b));
                return result;
            }

            inline scalar dd2V_dT2_p(const scalar p, const scalar T,
                                 const scalar z) const {
                const scalar V = z * RR * T / p;
                const scalar dV= dVdT_p(p, T, z);
                const scalar first =
                    (dV*(2*RR*V + 2*RR*b- da_dT) - dd2a_dT2*V + dd2a_dT2*b)
                    / (3*p *sqr(V) + (2*b*p - 2*RR*T)*V + a - 2*b*RR*T - 3*sqr(b)*p);
                const scalar second = (RR*sqr(V) + (2*b*RR-da_dT)*V + da_dT*b - RR*sqr(b))
                    * (dV*(6*p*V+2*b*p-2*RR*T) - 2*RR*V + da_dT - 2*RR*b)
                    / sqr(3*p*sqr(V) + (2*b*p-2*RR*T)*V + a - 2*b*RR*T - 3*sqr(b)*p);
                return first+second;
            }

            inline scalar dhdp_T(scalar p, scalar T, scalar W) const {
                const scalar z = this->Z(p, T);
                const scalar V = z * RR * T / p;
                const scalar dHdp_T = V-T*dVdT_p(p, T, z);
                return dHdp_T/this->R(W);
            }

            inline scalar dedp_T(scalar p, scalar T, scalar W) const {
                const scalar z = this->Z(p, T);
                [[maybe_unused]] const scalar V = z * RR * T / p;
                const scalar dUdp_T = -T*dVdT_p(p, T, z) - p*dVdp_T(p, T, z);
                return dUdp_T/this->R(W);
            }
        };

    public:
        // Constructors
        inline PengRobinsonGas(
            const Specie&sp,
            const realFluidProperty& prop);

        //- Construct from components
        inline PengRobinsonGas(
                    const Specie &sp
                );

        //- Construct from name and dictionary
        PengRobinsonGas(const word &name, const dictionary &dict);

        //- Construct as named copy
        inline PengRobinsonGas(const word &name, const PengRobinsonGas &);

        //- Construct and return a clone
        inline autoPtr<PengRobinsonGas> clone() const;


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
        inline scalar omega() const noexcept {
            return this->property_.omega_;
        }
        inline bool preferGas() const noexcept {
            return this->property_.is_gas;
        }

        //- Return the instantiated type name
        static word typeName() {
            return "PengRobinsonGas<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties
        inline EOSCore core(scalar p, scalar T) const {
            const auto Tc = this->Tc();
            const auto Pc = this->Pc();
            const scalar omega=this->omega();
            // const scalar Tr = T/Tc;
            const scalar a_c = 0.45724* sqr(RR*Tc)/Pc;
            const scalar kappa = 0.37464+1.54226*omega - 0.26992*sqr(omega);
            const scalar sqrtAlpha = 1+kappa*(1-sqrt(T/Tc));
            const scalar alpha= sqr(sqrtAlpha);
            const scalar dalpha_dT=-kappa* sqrt(alpha/(T*Tc));
            const scalar d2alpha_dT2 = kappa *(alpha-T*dalpha_dT)/(2*T*sqrt(T*Tc*alpha));

            const scalar a = a_c*alpha;
            const scalar b =0.077780*RR*Tc/Pc;
            const scalar da_dT = a_c*dalpha_dT;
            const scalar d2a_dT2 = a_c*d2alpha_dT2;

            return {a,b,da_dT,d2a_dT2,this->preferGas()};
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

        //- Return dCv/dT_p contribution [J/(kg K^2)]
        inline scalar dCvdT(const scalar p, const scalar T) const {
            return this->core(p,T).dCvdT(p,T,this->W());
        }

        //- Return Cv contribution [J/(kg K)]
        inline scalar Cv(scalar p, scalar T) const{
            return this->core(p,T).Cv(p,T,this->W());
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
            return this->core(p,T).psi(p,T, this->W());
        }

        //- Return (Cp - Cv) [J/(kg K]
        inline scalar CpMCv(scalar p, scalar T) const{
            return this->core(p,T).CpMCv(p,T, this->W());
        }

        //- Return volumetric coefficient of thermal expansion [1/T]
        inline scalar alphav(const scalar p, const scalar T) const{
            return this->core(p,T).alphav(p,T, this->W());
        }

        inline scalar dhdp_T(const scalar p, const scalar T) const {
            return this->core(p, T).dhdp_T(p, T, this->W());
        }

        inline scalar dedp_T(const scalar p, const scalar T) const {
            return this->core(p, T).dedp_T(p, T, this->W());
        }

        // mixing

        class EOSMixer {
        protected:
            const label nSpecie;
            // row major
            scalarSquareMatrix kij;
            explicit EOSMixer(label nSp):nSpecie(nSp), kij(nSp) { }
        public:
            EOSMixer(const speciesTable &spTable, const dictionary&physicalPropDict):
                EOSMixer(spTable.size()) {
                const dictionary* mixingRuleDict=physicalPropDict.subDictPtr("mixingRule");
                if(mixingRuleDict==nullptr) { // keep default, all k_ij=0
                    return;
                }
                const dictionary* binaryInterfactionDict=mixingRuleDict->subDictPtr("binaryInteraction");
                if(binaryInterfactionDict==nullptr) {
                    return;
                }
                const dictionary*kDict=binaryInterfactionDict->subDictPtr("k");
                if(kDict==nullptr) {
                    return;
                }
                const bool symmetry=kDict->lookupOrDefault<bool>("symmetry", true);

                auto skipFun=[](const word&key) {
                    if(key=="symmetry") {
                        return true;
                    }
                    return false;
                };
                parseBinaryInteractionMatrix(*kDict,skipFun,this->kij,spTable,symmetry);

                Info<<"Binary interaction coefficient table for PR equationOfState: \nk="<<this->kij<<endl;
            }

            inline scalar &k(label i, label j)& {
                return this->kij(i,j);
            }
            inline scalar k(label i,label j) const &{
                return this->kij(i,j);
            }

            template<class getCoreFun, class getXiFun>
            inline EOSCore mix(label nSpecie,
                               getCoreFun &coreFun,
                               getXiFun &specieMoleFracFun)const {

                scalar aSum=0;
                scalar da_dT_sum=0;
                scalar d2a_dT2_sum=0;
                scalar bSum=0;
                const bool is_gas=coreFun(0).is_gas;
                for(label i=0;i<nSpecie;i++) {
                    const scalar Xi=specieMoleFracFun(i);
                    const EOSCore iCore=coreFun(i);

                    bSum+=Xi*iCore.b;

                    for(label j=0;j<nSpecie;j++) {
                        const scalar Xj=specieMoleFracFun(j);
                        const EOSCore jCore=coreFun(j);
                        const scalar kij=this->k(i,j);
                        const scalar ai=iCore.a;
                        const scalar aj=jCore.a;
                        const scalar sqrt_aiaj=sqrt(ai * aj);
                        const scalar ddT_aiaj = ai*jCore.da_dT + aj*iCore.da_dT;
                        aSum += Xi*Xj*(1-kij)*sqrt_aiaj ;

                        // compute da/dT
                        const scalar ddT_sqrt_aiaj=1/(2*sqrt_aiaj)*ddT_aiaj;
                        da_dT_sum+= Xi*Xj*(1-kij) * ddT_sqrt_aiaj;
                        
                        // compute d2a/dT2
                        const scalar d2_dT2_sqrt_aiaj = 
                            1/(4*ai*aj * sqrt_aiaj) * (
                                2*sqr(ai)*aj *jCore.dd2a_dT2
                                - sqr(ai*jCore.da_dT)
                                + 2*ai *iCore.da_dT*aj*jCore.da_dT
                                + sqr(aj) * (2*ai *iCore.dd2a_dT2 - sqr(iCore.da_dT))
                            );
                        d2a_dT2_sum += Xi*Xj*(1-kij) *d2_dT2_sqrt_aiaj;
                    }
                }
                return EOSCore{aSum,bSum,da_dT_sum,d2a_dT2_sum,is_gas};
            }

        };

        // IO

        //- Write to Ostream
        void write(Ostream &os) const;

        // Member Operators

        inline void operator+=(const PengRobinsonGas &);

        inline void operator*=(const scalar);


        // Friend operators

        friend PengRobinsonGas operator+<Specie>
                (
                        const PengRobinsonGas &,
                        const PengRobinsonGas &
                );

        friend PengRobinsonGas operator*<Specie>
                (
                        const scalar s,
                        const PengRobinsonGas &
                );

        friend PengRobinsonGas operator==<Specie>
                (
                        const PengRobinsonGas &,
                        const PengRobinsonGas &
                );


        // Ostream Operator

        friend Ostream &operator<<<Specie>
                (
                        Ostream &,
                        const PengRobinsonGas &
                );
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PengRobinsonGasI.H"

#ifdef NoRepository
#include "PengRobinsonGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
