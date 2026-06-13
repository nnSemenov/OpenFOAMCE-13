/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "AttouFerschneider.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(AttouFerschneider, 0);
    addToRunTimeSelectionTable(dragModel, AttouFerschneider, dictionary);
}
}


Foam::List<Foam::word> read_solids(const Foam::dictionary& dict) {
    using Foam::word,Foam::List;

    if(dict.found("solid")) {
        word single=dict.lookup("solid");
        List ret{single};
        return ret;
    }

    return dict.lookup<List<word>>("solids");
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<const Foam::phaseModel*>
Foam::dragModels::AttouFerschneider::solidPhases(const phaseSystem&ps) const {
    Foam::List<const Foam::phaseModel*> ret(this->solidNames_.size());
    forAll(this->solidNames_, index) {
        ret[index]=&ps.phases()[this->solidNames_[index]];
    }
    return ret;
}

Foam::averageDiameter::sauterAverage Foam::dragModels::AttouFerschneider::
    get_solid_phase_info(const Foam::phaseSystem &ps) const
{
    List<const phaseModel*> solids= solidPhases(ps);
    auto get_phase_fun = [&](label i) -> const phaseModel & {
        return *solids[i];
    };

    return averageDiameter::sauter(solids.size(), get_phase_fun);
}

Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::KGasLiquid
(
    const phaseModel& gas,
    const phaseModel& liquid
) const
{
    averageDiameter::sauterAverage solid_info =
        get_solid_phase_info(gas.fluid());

    const volScalarField oneMinusGas(max(1 - gas, liquid.residualAlpha()));
    const volScalarField cbrtR
    (
        cbrt(solid_info.alpha/oneMinusGas)
    );
    const volScalarField magURel(mag(gas.U() - liquid.U()));

    return E1_ * gas.fluidThermo().mu() * sqr(oneMinusGas / solid_info.d) *
        sqr(cbrtR) / max(gas, gas.residualAlpha()) +
        E2_ * gas.rho() * magURel * (1 - gas) / solid_info.d * cbrtR;
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::KGasSolid
(
    const phaseModel& gas,
    const phaseModel& solid
) const
{
    const volScalarField oneMinusGas(max(1 - gas, solid.residualAlpha()));
    const volScalarField cbrtR
    (
        cbrt(max(solid, solid.residualAlpha())/oneMinusGas)
    );

    return
        E1_*gas.fluidThermo().mu()*sqr(oneMinusGas/solid.d())*sqr(cbrtR)
       /max(gas, gas.residualAlpha())
      + E2_*gas.rho()*mag(gas.U())*(1 - gas)/solid.d()*cbrtR;
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::KLiquidSolid
(
    const phaseModel& liquid,
    const phaseModel& solid
) const
{
    const phaseModel& gas = liquid.fluid().phases()[gasName_];

    return
        E1_*liquid.fluidThermo().mu()
       *sqr(max(solid, solid.residualAlpha())/solid.d())
       /max(liquid, liquid.residualAlpha())
      + E2_*liquid.rho()*mag(gas.U())*solid/solid.d();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::AttouFerschneider::AttouFerschneider
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dragModel(dict, interface, registerObject),
    interface_(interface),
    gasName_(dict.lookup("gas")),
    liquidName_(dict.lookup("liquid")),
    solidNames_(read_solids(dict)),
    E1_("E1", dimless, dict),
    E2_("E2", dimless, dict)
{
    if(solidNames_.empty()) {
        FatalErrorInFunction<<"Required at least 1 solid phase"<<exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::AttouFerschneider::~AttouFerschneider()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::K() const
{
    const phaseModel& gas = interface_.fluid().phases()[gasName_];
    const phaseModel& liquid = interface_.fluid().phases()[liquidName_];
    List<const phaseModel*> solids= solidPhases(interface_.fluid());

    const phaseModel* matched_solid_phase= nullptr;
    for(const phaseModel* solid: solids) {
        if(interface_.contains(*solid)) {
            matched_solid_phase=solid;
            break;
        }
    }
//    const phaseModel& solid = interface_.fluid().phases()[solidName_];

    if (interface_.contains(gas) && interface_.contains(liquid))
    {
        return KGasLiquid(gas, liquid);
    }
    if (interface_.contains(gas) && matched_solid_phase)
    {
        return KGasSolid(gas, *matched_solid_phase);
    }
    if (interface_.contains(liquid) && matched_solid_phase)
    {
        return KLiquidSolid(liquid, *matched_solid_phase);
    }

    FatalErrorInFunction
        << "The interface " << interface_.name() << " does not contain two "
        << "out of the gas, liquid and solid phase models."
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::AttouFerschneider::Kf() const
{
    return fvc::interpolate(K());
}


// ************************************************************************* //
