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
\*---------------------------------------------------------------------------*/

#include "idealLiquidMulticomponentMixture.H"

#include "../../specie/equationOfState/cubicEOS.H"

template <class ThermoType>
template <class Method, class... Args>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::massWeighted(Method psiMethod,
                                                 const Args &...args) const
{
    scalar psi = 0;

    forAll(Y_, i) { psi += Y_[i] * (specieThermos_[i].*psiMethod)(args...); }

    return psi;
}

template <class ThermoType>
template <class Method, class... Args>
Foam::scalar Foam::idealLiquidMulticomponentMixture<ThermoType>::
    thermoMixtureType::harmonicMassWeighted(Method psiMethod,
                                            const Args &...args) const
{
    scalar rPsi = 0;

    forAll(Y_, i) { rPsi += Y_[i] / (specieThermos_[i].*psiMethod)(args...); }

    return 1 / rPsi;
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::limit(const scalar T) const
{
    return T;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoType>
Foam::idealLiquidMulticomponentMixture<
    ThermoType>::idealLiquidMulticomponentMixture(const dictionary &dict)
    : multicomponentMixture<ThermoType>(dict),
      thermoMixture_(this->specieThermos()),
      transportMixture_(this->specieThermos())
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::W() const
{
    return harmonicMassWeighted(&ThermoType::W);
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::rho(scalar p, scalar T) const
{
    return harmonicMassWeighted(&ThermoType::rho, p, T);
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::psi(scalar p, scalar T) const
{
    return this->rho(p, T) / p;
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::alphav(scalar p, scalar T) const
{
    scalar volume_frac = 0, alphav_volume_frac = 0;
    forAll(Y_, i)
    {
        const scalar rhoi = specieThermos_[i].rho(p, T);
        const scalar Yi_over_rhoi = Y_[i] / rhoi;
        volume_frac += Yi_over_rhoi;
        alphav_volume_frac += specieThermos_[i].alphav(p, T) * Yi_over_rhoi;
    }

    return alphav_volume_frac / volume_frac;
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::hf() const
{
    return massWeighted(&ThermoType::hf);
}

#define thermoMixtureFunction(Func)                                            \
                                                                               \
    template <class ThermoType>                                                \
    Foam::scalar Foam::idealLiquidMulticomponentMixture<                       \
        ThermoType>::thermoMixtureType::Func(scalar p, scalar T) const         \
    {                                                                          \
        return massWeighted(&ThermoType::Func, p, T);                          \
    }

thermoMixtureFunction(Cp) thermoMixtureFunction(Cv) thermoMixtureFunction(hs)
    thermoMixtureFunction(ha) thermoMixtureFunction(es)
        thermoMixtureFunction(ea) thermoMixtureFunction(dhdp_T)
            thermoMixtureFunction(dedp_T)

                template <class ThermoType>
                Foam::scalar Foam::idealLiquidMulticomponentMixture<
                    ThermoType>::thermoMixtureType::dhedp_T(scalar p,
                                                            scalar T) const
{

    if (ThermoType::enthalpy()) { // can be constexpr on C++17
        return this->dhdp_T(p, T);
    }
    else {
        return this->dedp_T(p, T);
    }
}

// Cp or Cv
// thermoMixtureFunction(Cpv)

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::Cpv(scalar p, scalar T) const
{
    if (ThermoType::enthalpy()) { // can be constexpr on C++17
        return this->Cp(p, T);
    }
    else {
        return this->Cv(p, T);
    }
}
// hs, ha, es or ea
template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::he(scalar p, scalar T) const
{
    constexpr bool absolute = ThermoType::absolute();
    if (ThermoType::enthalpy()) { // can be constexpr on C++17
        if (absolute) {
            return this->ha(p, T);
        }
        else {
            return this->hs(p, T);
        }
    }
    else {
        if (absolute) {
            return this->ea(p, T);
        }
        else {
            return this->es(p, T);
        }
    }
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::gamma(scalar p, scalar T) const
{
    return this->Cp(p, T) / this->Cv(p, T);
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType::The(const scalar he, scalar p,
                                        scalar T0) const
{
    return ThermoType::T(*this, he, p, T0, &thermoMixtureType::he,
                         &thermoMixtureType::Cpv, &thermoMixtureType::limit);
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::transportMixtureType::mu(scalar p, scalar T) const
{
    auto prop = [](const ThermoType &specie, scalar p, scalar T) {
        const scalar mu = specie.mu(p, T);
        assert(mu > 0);
        return std::log(mu);
    };
    const scalar log_sum = moleWeighed(prop, p, T);
    return std::exp(log_sum);
}

template <class ThermoType>
Foam::scalar Foam::idealLiquidMulticomponentMixture<
    ThermoType>::transportMixtureType::kappa(scalar p, scalar T) const
{
    auto prop = [](const ThermoType &specie, scalar p, scalar T) {
        const scalar kappa = specie.kappa(p, T);
        assert(kappa > 0);
        return 1 / sqr(kappa);
    };
    const scalar kappa_to_neg2 = massWeighed(prop, p, T);
    return 1 / sqrt(kappa_to_neg2);
}

template <class ThermoType>
const typename Foam::idealLiquidMulticomponentMixture<
    ThermoType>::thermoMixtureType &
Foam::idealLiquidMulticomponentMixture<ThermoType>::thermoMixture(
    const scalarFieldListSlice &Y) const
{
    scalar sumX = 0;
    forAll(Y, i)
    {
        thermoMixture_.Y_[i] = Y[i];
        thermoMixture_.X_[i] = Y[i] / this->specieThermos()[i].W();
        sumX += thermoMixture_.X_[i];
    }

    forAll(Y, i) { thermoMixture_.X_[i] /= sumX; }

    return thermoMixture_;
}

template <class ThermoType>
const typename Foam::idealLiquidMulticomponentMixture<
    ThermoType>::transportMixtureType &
Foam::idealLiquidMulticomponentMixture<ThermoType>::transportMixture(
    const scalarFieldListSlice &Y) const
{
    scalar sumX = 0;

    forAll(Y, i)
    {
        transportMixture_.Y_[i] = Y[i];
        transportMixture_.X_[i] = Y[i] / this->specieThermos()[i].W();
        sumX += transportMixture_.X_[i];
    }

    forAll(Y, i) { transportMixture_.X_[i] /= sumX; }

    return transportMixture_;
}

template <class ThermoType>
const typename Foam::idealLiquidMulticomponentMixture<
    ThermoType>::transportMixtureType &
Foam::idealLiquidMulticomponentMixture<ThermoType>::transportMixture(
    const scalarFieldListSlice &Y, const thermoMixtureType &) const
{
    return transportMixture(Y);
}
