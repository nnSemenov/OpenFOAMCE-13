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
\*---------------------------------------------------------------------------*/

#include "porousMediaFluid.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include <cassert>
#include <optional>

using namespace Foam;

tmp<volScalarField::Internal> solvers::porousMediaFluid::
    TEqnNonIdealityPressureTerm(const volScalarField &dhedp_T) const
{
    auto left = rho_ * U_ & fvc::grad(p());
    return dhedp_T.internalField() *
        (left->internalField() + rho_.internalField() * dpdt);
}

fvScalarMatrix solvers::porousMediaFluid::TEqnCore(
    porousPhaseInfo::heatTransferInfoNode &node)
{
    auto getAlpha = [&](const word &phaseName) -> const volScalarField & {
        if (phaseName == porousPhaseInfo::FLUID_PHASE_NAME) {
            return this->porosity_;
        }
        return this->porousPhases.porousPhases().at(phaseName).alpha;
    };
    auto getRho = [&](const word &phaseName) -> const volScalarField & {
        if (phaseName == porousPhaseInfo::FLUID_PHASE_NAME) {
            return this->rho_;
        }
        return this->porousPhases.porousPhases().at(phaseName).thermo->rho();
    };
    auto getCpv = [&](const word &phaseName) -> tmp<volScalarField> {
        if (phaseName == porousPhaseInfo::FLUID_PHASE_NAME) {
            return this->thermoPtr_->Cpv();
        }
        return this->porousPhases.porousPhases().at(phaseName).thermo->Cpv();
    };
    auto getKappa = [&](const word &phaseName) -> tmp<volScalarField> {
        if (phaseName == porousPhaseInfo::FLUID_PHASE_NAME) {
            return this->thermophysicalTransport->kappaEff();
        }
        return this->porousPhases.porousPhases().at(phaseName).thermo->kappa();
    };

    volScalarField alphaRhoCp(
        IOobject("alphaRhoCpSum", this->runTime.name(), this->mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE,
                 false // no register
                 ),
        this->mesh_,
        dimensionedScalar("rhoCpEffTemp",
                          dimEnergy / dimVolume / dimTemperature, scalar(0)));
    volScalarField alphaSum(
        IOobject("alphaSum", this->runTime.name(), this->mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE,
                 false // no register
                 ),
        this->mesh_, dimensionedScalar("alphaSumTemp", dimless, scalar(0)));
    std::vector<effectiveHeatConductivity::kappaInfo> kappaInfos;
    for (const word &phaseName : node.thermalEquilibriumPhaseNames) {
        alphaRhoCp +=
            getAlpha(phaseName) * getRho(phaseName) * getCpv(phaseName);
        kappaInfos.emplace_back(getAlpha(phaseName), getKappa(phaseName));
        alphaSum += getAlpha(phaseName);
    }
    alphaSum.max(this->alpha_small_);
    tmp<volScalarField> kappaEff{nullptr};
    const bool singlePhase = node.thermalEquilibriumPhaseNames.size() == 1;
    if (singlePhase) {
        kappaEff = kappaInfos[0].kappa().clone();
    }
    else {
        kappaEff = node.effectiveKappaModel->kappaEff(kappaInfos);
    }
    assert(min(kappaEff.ref()).value() > 0);
    kappaEff.ref() *= alphaSum;
    assert(min(kappaEff.ref()).value() > 0);
    // Build fvMatrix
    volScalarField &T = getTFieldRef(node.thermalEquilibriumPhaseNames[0]);
    fvScalarMatrix TEqn(alphaRhoCp * fvm::ddt(T) - fvm::laplacian(kappaEff, T) -
                        fvModels().source(alphaRhoCp, T));
    // add fluid terms
    const std::optional<size_t> fluidIndex = node.fluidPhaseIndex();
    if (fluidIndex) {
        const volScalarField &he = thermoPtr_->he();
        const bool internalEnergy = he.name() == "e";
        const volScalarField &dhedp_T = thermoPtr_->dhedp_T();
        auto Cpv = this->thermoPtr_->Cpv();
        const surfaceScalarField phiCpv = (fvc::interpolate(Cpv) * phi_)();
        const auto &alphaF = this->porosity_;
        TEqn += fvm::div(phiCpv, T) - fvm::SuSp(fvc::div(phiCpv), T) +
            TEqnNonIdealityPressureTerm(dhedp_T);
        TEqn += alphaF * fvc::ddt(rho_, K) + fvc::div(phi_, K);
        TEqn += alphaF *
            pressureWork(internalEnergy ? fvc::div(phi_, p_() / rho_)()
                                        : -dpdt);

        if (buoyancy.valid()) {
            TEqn -= rho_ * (U_ & buoyancy->g);
        }
    }
    // add phasewise heat transfer
    for (auto &htSource : node.heatTransferSources) {
        const word phaseDest =
            node.thermalEquilibriumPhaseNames[htSource.phaseA_index];
        const word phaseSrc = htSource.phaseB_name;
        const volScalarField &Tdest = getTFieldRef(phaseDest);
        const volScalarField &Tsrc = getTFieldRef(phaseSrc);
        const auto &htInfo = htSource.heatTransferInfo;
        const auto hsAs = getAlpha(phaseSrc) * getAlpha(phaseDest) *
            htInfo.heatTransferCoefficient *
            htInfo.effectiveSpecificSurfaceArea;
        if (htInfo.explicitTerm) {
            TEqn -= hsAs * (Tsrc - Tdest);
        }
        else {
            TEqn -= hsAs.ref() * Tsrc;
            TEqn += fvm::Sp(hsAs, Tdest);
        }
    }

    return TEqn;
}

void solvers::porousMediaFluid::thermophysicalPredictor()
{
    auto htGraph = this->porousPhases.heatTransferSummary();
    const size_t N_eqn = htGraph.size();
    std::vector<fvScalarMatrix> TEqnList;
    for (auto &node : htGraph) {
        TEqnList.emplace_back(this->TEqnCore(node));
    }
    assert(TEqnList.size() == N_eqn);
    for (size_t i = 0; i < N_eqn; i++) {
        fvScalarMatrix &TEqn = TEqnList[i];
        TEqn.relax();
        fvConstraints().constrain(TEqn);
        TEqn.solve();

        volScalarField &T =
            getTFieldRef(htGraph[i].thermalEquilibriumPhaseNames[0]);
        fvConstraints().constrain(T);
        const auto &p_ = this->p_();
        for (size_t phIdx = 0;
             phIdx < htGraph[i].thermalEquilibriumPhaseNames.size(); phIdx++) {
            const word &phaseName =
                htGraph[i].thermalEquilibriumPhaseNames[phIdx];

            basicThermo &thermo = getThermo(phaseName);
            if (phIdx > 0) {
                thermo.T() = T;
            }
            thermo.he() = thermo.he(p_, T);
            thermo.correct();
        }
    }

    // auto & thermo=this->thermoPtr_();
    // volScalarField & T=thermo.T();
    // const volScalarField& he = thermo.he();
    // const bool internalEnergy = he.name() == "e";
    // const auto Cpv=thermo.Cpv();
    // const surfaceScalarField phiCpv =(fvc::interpolate(Cpv) * phi_)();

    // const volScalarField & dhedp_T = thermo.dhedp_T();

    // const volScalarField & alphaF = this->porosity_;

    // volScalarField kappaEff = (this->thermophysicalTransport->kappaEff() *
    // alphaF)(); volScalarField rhoCpEff = (alphaF* rho_ * Cpv)(); for(auto &
    // pair: this->porousPhases.porousPhases()) {
    //     const auto & alphaP=pair.second.alpha;
    //     auto & thermoP=pair.second.thermo;
    //     const volScalarField& rhoP=thermoP->rho();
    //     const volScalarField& CpvP=thermoP->Cpv();
    //     rhoCpEff+=alphaP*rhoP*CpvP;
    //     kappaEff+=alphaP*thermoP->kappa();
    // }
    // // Thermal equilbrium
    // fvScalarMatrix TEqn (

    //     rhoCpEff* fvm::ddt(T)
    //     + alphaF*(fvm::div(phiCpv, T)
    //         - fvm::SuSp(fvc::div(phiCpv), T))
    //     + alphaF*TEqnNonIdealityPressureTerm(dhedp_T)

    //     + alphaF*(fvc::ddt(rho_, K) + fvc::div(phi_, K)) // Fluid kinetic
    //     - fvm::laplacian(kappaEff, T)
    //     + alphaF*pressureWork
    //       (
    //         internalEnergy
    //         ? fvc::div(phi_, p_()/rho_)()
    //         : -dpdt
    //       )
    //     ==
    //     (
    //         buoyancy.valid()
    //     ? fvModels().source(rhoCpEff, T) + alphaF*rho_*(U_ & buoyancy->g)
    //     : fvModels().source(rhoCpEff, T)
    //     )
    // );

    // Fluid only
    // fvScalarMatrix TEqn (
    //     rho_ * Cpv* fvm::ddt(T)
    //     + fvm::div(phiCpv, T)
    //         - fvm::SuSp(fvc::div(phiCpv), T)
    //     + TEqnNonIdealityPressureTerm(dhedp_T)

    //     + fvc::ddt(rho_, K) + fvc::div(phi_, K)
    //     - fvm::laplacian(kappaEff, T)
    //     + pressureWork
    //       (
    //         internalEnergy
    //         ? fvc::div(phi_, p_()/rho_)()
    //         : -dpdt
    //       )
    //     ==
    //         (
    //             buoyancy.valid()
    //         ? fvModels().source(rho_*Cpv, T) + rho_*(U_ & buoyancy->g)
    //         : fvModels().source(rho_*Cpv, T)
    //         )
    // );

    // TEqn.relax();

    // fvConstraints().constrain(TEqn);

    // TEqn.solve();

    // fvConstraints().constrain(T);
    // // Update energy
    // thermo.he() = thermo.he(this->p_(), T);
    // // Update other properties
    // this->thermoPtr_().correct();

    // for(auto & pair: this->porousPhases.porousPhases()) {
    //     auto & thermoP=pair.second.thermo;
    //     thermoP->T()=T;
    //     thermoP->he() = thermoP->he(this->p_(), T);
    //     thermoP->correct();
    // }
}