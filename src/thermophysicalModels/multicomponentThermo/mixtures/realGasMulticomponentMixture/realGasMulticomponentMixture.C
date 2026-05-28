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

#include "realGasMulticomponentMixture.H"
#include "../../specie/equationOfState/cubicEOS.H"

template <class ThermoType>
template <class Method, class... Args>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::massWeighted(Method psiMethod,
                                                 const Args &...args) const
{
    scalar psi = 0;

    forAll(Y_, i) { psi += Y_[i] * (specieThermos_[i].*psiMethod)(args...); }

    return psi;
}

template <class ThermoType>
template <class Method, class... Args>
Foam::scalar Foam::realGasMulticomponentMixture<ThermoType>::thermoMixtureType::
    harmonicMassWeighted(Method psiMethod, const Args &...args) const
{
    scalar rPsi = 0;

    forAll(Y_, i) { rPsi += Y_[i] / (specieThermos_[i].*psiMethod)(args...); }

    return 1 / rPsi;
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::limit(const scalar T) const
{
    return T;
}

template <class ThermoType>
template <class Method, class... Args>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::transportMixtureType::moleWeighted(Method psiMethod,
                                                    const Args &...args) const
{
    scalar psi = 0;

    forAll(X_, i) { psi += X_[i] * (specieThermos_[i].*psiMethod)(args...); }

    return psi;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoType>
Foam::realGasMulticomponentMixture<ThermoType>::realGasMulticomponentMixture(
    const dictionary &dict)
    : multicomponentMixture<ThermoType>(dict),
      thermoMixture_(this->specieThermos()),
      transportMixture_(this->specieThermos())
{
    static_assert(Foam::is_cubic_EOS<ThermoType>::value,
                  "Only cubic EOS is supported now.");

    this->thermoMixture_.mixer_ =
        new typename ThermoType::EOSMixer{this->species(), dict};
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class ThermoType>
// Foam::cubicEOSRawCoefficient
// Foam::realGasMulticomponentMixture<ThermoType>::thermoMixtureType::cubicEOSCoeffOfSpecie(scalar
// p, scalar T,
//                                                                                          label idx) const {
//     return
// }

// template<class ThermoType>
// Foam::cubicEOSRawCoefficient
// Foam::realGasMulticomponentMixture<ThermoType>::thermoMixtureType::mixedRawCoeff(scalar
// p, scalar T) const {
//
//     return
// }

template <class ThermoType>
auto Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::mixedCore(scalar p, scalar T) const
{

    auto coeffFun = [this, p, T](label idx) {
        return this->specieThermos_[idx].core(p, T);
    };
    auto XiFun = [this](label idx) { return this->X_[idx]; };

    return this->mixer_().mix(this->X_.size(), coeffFun, XiFun);
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::W() const
{
    return harmonicMassWeighted(&ThermoType::W);
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::rho(scalar p, scalar T) const
{
    auto mixedCore = this->mixedCore(p, T);
    return mixedCore.rho(p, T, this->W());
    //    return harmonicMassWeighted(&ThermoType::rho, p, T);
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::psi(scalar p, scalar T) const
{
    auto mixedCore = this->mixedCore(p, T);
    return mixedCore.psi(p, T, this->W());
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::hf() const
{
    return massWeighted(&ThermoType::hf);
}

#define thermoMixtureFunction(Func)                                            \
                                                                               \
    template <class ThermoType>                                                \
    Foam::scalar                                                               \
    Foam::realGasMulticomponentMixture<ThermoType>::thermoMixtureType::Func(   \
        scalar p, scalar T) const                                              \
    {                                                                          \
        return massWeighted(&ThermoType::Func, p, T);                          \
    }

#define realThermoMixtureFunction(Func, eosFunc)                               \
    template <class ThermoType>                                                \
    Foam::scalar                                                               \
    Foam::realGasMulticomponentMixture<ThermoType>::thermoMixtureType::Func(   \
        scalar p, scalar T) const                                              \
    {                                                                          \
        const scalar ideal = massWeighted(&ThermoType::ideal_##Func, p, T);    \
        const auto core = this->mixedCore(p, T);                               \
        const scalar residual = core.eosFunc(p, T, this->W());                 \
        return ideal + residual;                                               \
    }

// thermoMixtureFunction(Cp)
// thermoMixtureFunction(Cv)
// thermoMixtureFunction(hs)
// thermoMixtureFunction(ha)

realThermoMixtureFunction(Cp, Cp) realThermoMixtureFunction(Cv, Cv)
    realThermoMixtureFunction(hs, h) realThermoMixtureFunction(ha, h)
        realThermoMixtureFunction(es, e) realThermoMixtureFunction(ea, e)

            template <class ThermoType>
            Foam::scalar
    Foam::realGasMulticomponentMixture<ThermoType>::thermoMixtureType::dhdp_T(
        scalar p, scalar T) const
{
    const auto core = this->mixedCore(p, T);
    return core.dhdp_T(p, T, this->W());
}
template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::dedp_T(scalar p, scalar T) const
{
    const auto core = this->mixedCore(p, T);
    return core.dedp_T(p, T, this->W());
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::dhedp_T(scalar p, scalar T) const
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
Foam::scalar Foam::realGasMulticomponentMixture<
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
Foam::scalar Foam::realGasMulticomponentMixture<
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
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::gamma(scalar p, scalar T) const
{
    return this->Cp(p, T) / this->Cv(p, T);
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType::The(const scalar he, scalar p,
                                        scalar T0) const
{
    return ThermoType::T(*this, he, p, T0, &thermoMixtureType::he,
                         &thermoMixtureType::Cpv, &thermoMixtureType::limit);
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::transportMixtureType::mu(scalar p, scalar T) const
{
    return moleWeighted(&ThermoType::mu, p, T);
}

template <class ThermoType>
Foam::scalar Foam::realGasMulticomponentMixture<
    ThermoType>::transportMixtureType::kappa(scalar p, scalar T) const
{
    return moleWeighted(&ThermoType::kappa, p, T);
}

template <class ThermoType>
const typename Foam::realGasMulticomponentMixture<
    ThermoType>::thermoMixtureType &
Foam::realGasMulticomponentMixture<ThermoType>::thermoMixture(
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
const typename Foam::realGasMulticomponentMixture<
    ThermoType>::transportMixtureType &
Foam::realGasMulticomponentMixture<ThermoType>::transportMixture(
    const scalarFieldListSlice &Y) const
{
    scalar sumX = 0;

    forAll(Y, i)
    {
        transportMixture_.X_[i] = Y[i] / this->specieThermos()[i].W();
        sumX += transportMixture_.X_[i];
    }

    forAll(Y, i) { transportMixture_.X_[i] /= sumX; }

    return transportMixture_;
}

template <class ThermoType>
const typename Foam::realGasMulticomponentMixture<
    ThermoType>::transportMixtureType &
Foam::realGasMulticomponentMixture<ThermoType>::transportMixture(
    const scalarFieldListSlice &Y, const thermoMixtureType &) const
{
    return transportMixture(Y);
}
