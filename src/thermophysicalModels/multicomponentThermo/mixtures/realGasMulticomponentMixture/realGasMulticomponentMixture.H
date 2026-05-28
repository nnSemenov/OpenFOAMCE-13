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

Class
    Foam::realGasMulticomponentMixture

Description
    Thermophysical properties mixing class for real gas mixtures.
    Mixed equation of state are used for density, enthalpy, specific heat, etc.
    Thermophysical properties are calculated in rigorous way considering
    residual properties. Currently only cubic EOSs are supported.

    Mole-weighed are used for kappa and mu.
    TODO: Use better/various mixing rule for transport features.

SourceFiles
    realGasMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef realGasMulticomponentMixture_H
#define realGasMulticomponentMixture_H

#include "multicomponentMixture.H"
#include "FieldListSlice.H"

namespace Foam
{

// forward declaration
struct cubicEOSRawCoefficient;

template <class ThermoType>
class realGasMulticomponentMixture : public multicomponentMixture<ThermoType>
{
  public:
    // Public Classes

    //- Mixing type for thermodynamic properties
    class thermoMixtureType
    {
        // Private Data

        //- List of specie thermo
        const PtrList<ThermoType> &specieThermos_;

        //- List of mass fractions
        mutable List<scalar> Y_;
        //- List of mole fractions
        mutable List<scalar> X_;

        autoPtr<typename ThermoType::EOSMixer> mixer_;

        //- Calculate a mass-fraction-weighted property
        template <class Method, class... Args>
        scalar massWeighted(Method psiMethod, const Args &...args) const;

        //- Calculate a harmonic mass-fraction-weighted property
        template <class Method, class... Args>
        scalar harmonicMassWeighted(Method psiMethod,
                                    const Args &...args) const;

        //- Limit the given temperature
        scalar limit(const scalar T) const;

      protected:
        // cubicEOS properties
        //        cubicEOSRawCoefficient cubicEOSCoeffOfSpecie(scalar p, scalar
        //        T, label idx)const;
        //
        //        cubicEOSRawCoefficient mixedRawCoeff(scalar p, scalar T)
        //        const;

        auto mixedCore(scalar p, scalar T) const;

      public:
        friend class realGasMulticomponentMixture;

        // Constructors

        //- Construct from list of specie thermo
        thermoMixtureType(const PtrList<ThermoType> &specieThermos)
            : specieThermos_(specieThermos), Y_(specieThermos.size()),
              X_(specieThermos.size()), mixer_{nullptr}
        {}

        // Fundamental properties

        //- Molecular weight [kg/kmol]
        scalar W() const;

        //- Return density [kg/m^3]
        scalar rho(scalar p, scalar T) const;

        //- Return compressibility [s^2/m^2]
        scalar psi(scalar p, scalar T) const;

        // Heat capacity at constant pressure [J/kg/K]
        scalar Cp(const scalar p, const scalar T) const;

        // Heat capacity at constant volume [J/kg/K]
        scalar Cv(const scalar p, const scalar T) const;

        // Sensible enthalpy [J/kg]
        scalar hs(const scalar p, const scalar T) const;

        // Absolute enthalpy [J/kg]
        scalar ha(const scalar p, const scalar T) const;

        // Enthalpy of formation [J/kg]
        scalar hf() const;

        // Sensible internal energy [J/kg]
        scalar es(const scalar p, const scalar T) const;

        // Absolute internal energy [J/kg]
        scalar ea(const scalar p, const scalar T) const;

        // Mass specific derived properties

        //- Heat capacity at constant pressure/volume [J/kg/K]
        scalar Cpv(const scalar p, const scalar T) const;

        //- Gamma = Cp/Cv []
        scalar gamma(const scalar p, const scalar T) const;

        //- Enthalpy/Internal energy [J/kg]
        scalar he(const scalar p, const scalar T) const;

        //- Return isothermal enthalpy-pressure derivative [m^3/kg]
        scalar dhdp_T(const scalar p, const scalar T) const;

        //- Return isothermal internal energy-pressure derivative [m^3/kg]
        scalar dedp_T(const scalar p, const scalar T) const;

        //- Return isothermal enthalpy/internal-energy pressure derivative
        //[m^3/kg]
        scalar dhedp_T(const scalar p, const scalar T) const;

        // Energy->temperature  inversion functions

        //- Temperature from enthalpy or internal energy
        //  given an initial temperature T0
        scalar The(const scalar he, const scalar p, const scalar T0) const;
    };

    //- Mixing type for transport properties
    class transportMixtureType
    {
        // Private Data

        //- List of specie thermo
        const PtrList<ThermoType> &specieThermos_;

        //- List of mole fractions
        mutable List<scalar> X_;

        //- Calculate a mole-fraction-weighted property
        template <class Method, class... Args>
        scalar moleWeighted(Method psiMethod, const Args &...args) const;

      public:
        friend class realGasMulticomponentMixture;

        // Constructors

        //- Construct from list of specie thermo
        transportMixtureType(const PtrList<ThermoType> &specieThermos)
            : specieThermos_(specieThermos), X_(specieThermos.size())
        {}

        // Transport properties

        //- Dynamic viscosity [kg/m/s]
        scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        scalar kappa(const scalar p, const scalar T) const;
    };

  private:
    // Private Data

    //- Mutable storage for the cell/face mixture thermo data
    mutable thermoMixtureType thermoMixture_;

    //- Mutable storage for the cell/face mixture transport data
    mutable transportMixtureType transportMixture_;

  public:
    // Constructors

    //- Construct from a dictionary
    realGasMulticomponentMixture(const dictionary &);

    //- Disallow default bitwise copy construction
    realGasMulticomponentMixture(
        const realGasMulticomponentMixture<ThermoType> &) = delete;

    //- Destructor
    virtual ~realGasMulticomponentMixture() {}

    // Member Functions

    //- Return the instantiated type name
    static word typeName()
    {
        return "realGasMulticomponentMixture<" + ThermoType::typeName() + '>';
    }

    //- Return the mixture for thermodynamic properties
    const thermoMixtureType &thermoMixture(const scalarFieldListSlice &) const;

    //- Return the mixture for transport properties
    const transportMixtureType &transportMixture(
        const scalarFieldListSlice &) const;

    //- Return the mixture for transport properties
    const transportMixtureType &transportMixture(
        const scalarFieldListSlice &, const thermoMixtureType &) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "realGasMulticomponentMixture.C"
#endif

#endif // realGasMulticomponentMixture_H
