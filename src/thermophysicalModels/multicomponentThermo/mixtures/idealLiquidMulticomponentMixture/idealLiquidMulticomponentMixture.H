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
    Foam::idealLiquidMulticomponentMixture

Description
    Thermophysical properties mixing class for ideal liquid mixture.
    Apply harmonic mass weighed mixing for density (no extra volume),
    mass weighed mixing for he and Cp Cv (no extra energy),
    Arrhenius mixing for mu, Vredeveled mixing for kappa (n=-2).

    That is:
    log(mu) = sum n_i * log(mu_i)
    kappa^-2 = sum w_i * kappa_i^-2

SourceFiles
    idealLiquidMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef idealLiquidMulticomponentMixture_H
#define idealLiquidMulticomponentMixture_H

#include "multicomponentMixture.H"
#include "FieldListSlice.H"

namespace Foam
{

template <class ThermoType>
class idealLiquidMulticomponentMixture
    : public multicomponentMixture<ThermoType>
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

        //- Calculate a mass-fraction-weighted property
        template <class Method, class... Args>
        scalar massWeighted(Method psiMethod, const Args &...args) const;

        //- Calculate a harmonic mass-fraction-weighted property, also known as
        // volume fraction weighted
        template <class Method, class... Args>
        scalar harmonicMassWeighted(Method psiMethod,
                                    const Args &...args) const;

        //- Limit the given temperature
        scalar limit(const scalar T) const;

      protected:
      public:
        friend class idealLiquidMulticomponentMixture;

        // Constructors

        //- Construct from list of specie thermo
        thermoMixtureType(const PtrList<ThermoType> &specieThermos)
            : specieThermos_(specieThermos), Y_(specieThermos.size()),
              X_(specieThermos.size())
        {}

        // Fundamental properties

        //- Molecular weight [kg/kmol]
        scalar W() const;

        //- Return density [kg/m^3]
        scalar rho(scalar p, scalar T) const;

        //- Return compressibility [s^2/m^2]
        scalar psi(scalar p, scalar T) const;

        //- Return volumetric coefficient of thermal expansion [1/T]
        scalar alphav(const scalar p, const scalar T) const;

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
        //- List of mass fractions
        mutable List<scalar> Y_;

        //        mutable List<scalar> PHI_;

        template <class Method, class... Args>
        inline scalar moleWeighed(const Method &prop, const Args &...args) const
        {
            scalar sum = 0;
            forAll(this->X_, i)
            {
                sum += X_[i] * (prop(this->specieThermos_[i], args...));
            }
            return sum;
        }

        template <class Method, class... Args>
        inline scalar massWeighed(const Method &prop, const Args &...args) const
        {
            scalar sum = 0;
            forAll(this->Y_, i)
            {
                sum += Y_[i] * (prop(this->specieThermos_[i], args...));
            }
            return sum;
        }

      public:
        friend class idealLiquidMulticomponentMixture;

        // Constructors

        //- Construct from list of specie thermo
        transportMixtureType(const PtrList<ThermoType> &specieThermos)
            : specieThermos_(specieThermos), X_(specieThermos.size()),
              Y_(specieThermos.size())
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
    idealLiquidMulticomponentMixture(const dictionary &);

    //- Disallow default bitwise copy construction
    idealLiquidMulticomponentMixture(
        const idealLiquidMulticomponentMixture<ThermoType> &) = delete;

    //- Destructor
    virtual ~idealLiquidMulticomponentMixture() {}

    // Member Functions

    //- Return the instantiated type name
    static word typeName()
    {
        return "idealLiquidMulticomponentMixture<" + ThermoType::typeName() +
            '>';
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
#include "idealLiquidMulticomponentMixture.C"
#endif

#endif // idealLiquidMulticomponentMixture_H
