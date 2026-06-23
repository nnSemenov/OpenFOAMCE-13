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

#include <cstdlib>

#include "uniformDimensionedFields.H"
#include "Niegodajew.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam::heatTransferModels
{
defineTypeNameAndDebug(Niegodajew, 0);
addToRunTimeSelectionTable(heatTransferModel, Niegodajew, dictionary);
} // namespace Foam::heatTransferModels

Foam::heatTransferModels::Niegodajew::Niegodajew(
    const dictionary &dict, const phaseInterface &interface,
    const bool registerObject)
    : heatTransferModel{dict, interface, registerObject},
      interface_(
          interface.modelCast<heatTransferModel, dispersedPhaseInterface>()),
      gasName_{dict.lookup<word>("gas")},
      liquidName_{dict.lookup<word>("liquid")},
      solidNames_{dict.lookup<List<word>>("solids")},
      residualRe_{dict.lookupOrDefault<scalar>("residualRe", 1e-6)},
      sigma_{"sigma", dimForce / dimLength, dict}
{

    if (solidNames_.empty()) {
        FatalErrorInFunction << "Required at least 1 solid phase"
                             << exit(FatalError);
    }
    if (residualRe_ <= 0) {
        FatalErrorInFunction << "Residual Re must be positive, found "
                             << residualRe_ << exit(FatalError);
    }

    const word continuous_name = interface_.continuous().name();
    const word dispersed_name = interface_.dispersed().name();

    const bool phase_ok = [&]() {
        if (gasName_ == continuous_name and liquidName_ == dispersed_name) {
            return true;
        }
        if (gasName_ == dispersed_name and liquidName_ == continuous_name) {
            return true;
        }
        return false;
    }();
    if (not phase_ok) {
        FatalErrorInFunction << "Invalid phase model: gas = " << gasName_
                             << ", liquid = " << liquidName_
                             << ", continuous = " << continuous_name
                             << ", dispersed = " << dispersed_name
                             << ". Niegodajew is only valid for gas-liquid "
                                "interface heat transfer in trickle bed"
                             << exit(FatalError);
    }
}

Foam::averageDiameter::sauterAverage Foam::heatTransferModels::Niegodajew::
    average_diameter() const
{
    auto getPhaseModel = [&](label i) -> const phaseModel & {
        auto model_ptr =
            interface_.fluid().phases().lookup(this->solidNames_[i]);
        if (model_ptr == nullptr) {
            FatalErrorInFunction
                << "Solid phase " << this->solidNames_[i]
                << " is missing from phase model. Required by Niegodajew"
                << exit(FatalError);
            std::abort();
        }

        return *model_ptr;
    };

    return averageDiameter::sauter(this->solidNames_.size(), getPhaseModel);
}

Foam::tmp<Foam::volScalarField> Foam::heatTransferModels::Niegodajew::K(
    const scalar residualAlpha) const
{
    const auto [alpha_solid, dp] = average_diameter();

    assert(max(alpha_solid).value() > 0);
    assert(max(dp).value() > 0);

    const volScalarField porosity(1 - alpha_solid);
    assert(max(porosity).value() < 1);
    assert(min(porosity).value() > 0);

    // area/volume ratio
    const volScalarField a(6 * alpha_solid / dp);
    // effective diameter
    const volScalarField de(1. / a);
    assert(max(de).value() > 0);

    const phaseModel &liquid = interface_.fluid().phases()[liquidName_];
    const phaseModel &gas = interface_.fluid().phases()[gasName_];

    // Re_G = alpha_G * U_G * rho_G * dp / (1-porosity)/mu_G. The paper uses
    // superficial velocity for Re_G, alpha*U is equal.
    const volScalarField Re_G(max(static_cast<const volScalarField &>(gas) *
                                      mag(gas.U()) * gas.rho() * dp /
                                      (gas.fluidThermo().mu() * alpha_solid),
                                  residualRe_));
    assert(Re_G.dimensions() == dimless);
    assert(min(Re_G).value() > 0);

    const auto &g =
        interface_.mesh().lookupObject<uniformDimensionedVectorField>("g");

    const volScalarField Ga_G(sqr(gas.rho()) * mag(g) *
                              pow(de * porosity / alpha_solid, 3) /
                              sqr(gas.fluidThermo().mu()));
    assert(Ga_G.dimensions() == dimless);

    const volScalarField Eo(liquid.rho() * mag(g) *
                            sqr(de * porosity / alpha_solid) / sigma_);
    assert(Eo.dimensions() == dimless);
    // Note: this Nu is different from traditional OpenFOAM convention
    const volScalarField Nu(
        IOobject("Nu_" + interface_.name(), interface_.mesh(),
                 IOobject::NO_READ, IOobject::AUTO_WRITE, true),
        pow(Re_G, 1.169) * pow(Ga_G, -0.8399) * pow(Eo, 0.7176));

    if (interface_.mesh().time().writeTime()) {
        Nu.write();
    }

    // Definition of Nu differs from traditional correlations. Usually gas is
    // continuous phase, but here use liquid kappa and particle diameter. In
    // traditional correlations, K = 6 * alpha_disp * kappa_cont / d_disp^2
    // *Nu_trad; while Nu_trad = h * d_disp / kappa_cont.
    //
    // So, K = 6*h*alpha_disp/d_disp. To keep this convention, h = Nu_this *
    // kappa_liquid / d_e, so K =6 * alpha_disp/d_disp * Nu * k_liq / d_e

    return 6 * max(interface_.dispersed(), residualAlpha) * Nu *
        liquid.thermo().kappa() / (interface_.dispersed().d() * de);
}
