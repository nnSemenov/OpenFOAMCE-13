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