//
// Created by joseph on 25-7-11.
//

#ifndef RedlichKwongGas_H
#define RedlichKwongGas_H

#include "autoPtr.H"
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
    // private data
    //- Critical Temperature [K]
    scalar Tc_;

    //- Critical volume [m^3/kmol]
    scalar Vc_;

    //- Critical compression factor [-]
    scalar Zc_;

    //- Critical Pressure [Pa]
    scalar Pc_;

public:
    inline RedlichKwongGas(
            const Specie &sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Zc,
            const scalar& Pc
       );

    RedlichKwongGas(const word&name, const dictionary&dict);

    inline RedlichKwongGas(const word&name, const RedlichKwongGas&);

    inline autoPtr<RedlichKwongGas> clone() const;


    // Member Functions

    //- Return the instantiated type name
    static word typeName()
    {
        return "RedlichKwongGas<" + word(Specie::typeName_()) + '>';
    }

    // Fundamental properties

    inline Foam::cubicEOSCoefficient coeff_AB(scalar p, scalar T) const;

    //- Is the equation of state is incompressible i.e. rho != f(p)
    static const bool incompressible = false;

    //- Is the equation of state is isochoric i.e. rho = const
    static const bool isochoric = false;

    //- Return density [kg/m^3]
    inline scalar rho(scalar p, scalar T) const;

    //- Return enthalpy contribution [J/kg]
    inline scalar h(const scalar p, const scalar T) const;

    //- Return Cp contribution [J/(kg K]
    inline scalar Cp(scalar p, scalar T) const;

    //- Return internal energy contribution [J/kg]
    inline scalar e(const scalar p, const scalar T) const;

    //- Return Cv contribution [J/(kg K]
    inline scalar Cv(scalar p, scalar T) const;

    //- Return entropy contribution to the integral of Cp/T [J/kg/K]
    inline scalar sp(const scalar p, const scalar T) const;

    //- Return entropy contribution to the integral of Cv/T [J/kg/K]
    inline scalar sv(const scalar p, const scalar T) const;

    //- Return compressibility [s^2/m^2]
    inline scalar psi(scalar p, scalar T) const;

    //- Return compression factor []
    inline scalar Z(scalar p, scalar T) const;

    //- Return (Cp - Cv) [J/(kg K]
    inline scalar CpMCv(scalar p, scalar T) const;

    //- Return volumetric coefficient of thermal expansion [1/T]
    inline scalar alphav(const scalar p, const scalar T) const;


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
