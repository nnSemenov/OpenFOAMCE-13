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

Class
    Foam::solvers::porousMediaFluid

Description
    Solver module for compressible flow in porous media with heat transfer.
    Supports arbitrary number of porous phases, heat transfer between
    fluid-porous and porous-porous, both supporting thermal equilibrium and
    non-equilibrium.

SourceFiles
    porousMediaFluid.C
    thermophysicalPredictor.C
    porousPhaseInfo.H
    porousPhaseInfo.C
    effectiveHeatConductivity.H
    effectiveHeatConductivity.C

\*---------------------------------------------------------------------------*/

#ifndef porousMediaFluid_H
#define porousMediaFluid_H

#include "fluid.H"
#include "porousPhaseInfo.H"

namespace Foam
{

namespace solvers
{

class porousMediaFluid : public fluid
{
  protected:
    volVectorField U_physical;
    volScalarField porosity_;
    porousPhaseInfo porousPhases;

    scalar alpha_small_{1e-6};

    void updatePorosity();

    tmp<volScalarField::Internal> TEqnNonIdealityPressureTerm(
        const volScalarField &dhedp_T) const;

    fvScalarMatrix TEqnCore(porousPhaseInfo::heatTransferInfoNode &node);

    volScalarField &getTFieldRef(const word &phaseName)
    {
        if (phaseName == porousPhaseInfo::FLUID_PHASE_NAME) {
            return this->thermoPtr_->T();
        }
        return this->porousPhases.porousPhases().at(phaseName).thermo->T();
    }

    basicThermo &getThermo(const word &phaseName)
    {
        if (phaseName == porousPhaseInfo::FLUID_PHASE_NAME) {
            return this->thermoPtr_();
        }
        return this->porousPhases.porousPhases().at(phaseName).thermo();
    }

  public:
    TypeName("porousMediaFluid");

    // Constructors
    porousMediaFluid(fvMesh &mesh);

    porousMediaFluid(const porousMediaFluid &) = delete;

    // tmp<volScalarField> porosity() const;

    ~porousMediaFluid() override;

    void momentumPredictor() override;

    void thermophysicalPredictor() override;

    // void prePredictor() override;

    // void postCorrector() override;

    //- Disallow default bitwise assignment
    void operator=(const porousMediaFluid &) = delete;
};

} // namespace solvers
} // namespace Foam

#endif