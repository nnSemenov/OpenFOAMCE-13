#ifndef effectiveHeatConductivity_H
#define effectiveHeatConductivity_H

#include <span>
#include "volFieldsFwd.H"
#include "typeInfo.H"

namespace Foam
{

class effectiveHeatConductivity
{
  protected:
    scalar alpha_small_{1e-6};

  public:
    // TypeName("effectiveHeatConducitvity");

    explicit effectiveHeatConductivity(const dictionary &dict) {}

    virtual ~effectiveHeatConductivity() = default;

    struct kappaInfo
    {
        const volScalarField &alpha;
        tmp<volScalarField> kappa;
    };

    virtual tmp<volScalarField> kappaEff(std::span<kappaInfo>) const = 0;

    scalar alpha_small() const { return this->alpha_small_; }

    scalar &alpha_small() { return this->alpha_small_; }

    static autoPtr<effectiveHeatConductivity> New(const dictionary &dict);
};

} // namespace Foam

#endif // effectiveHeatConductivity_H
