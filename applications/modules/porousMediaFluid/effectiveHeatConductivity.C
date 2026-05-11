#include "effectiveHeatConductivity.H"
#include <cassert>
#include <cstdlib>

#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

using namespace Foam;

class parallel : public effectiveHeatConductivity
{
  public:
    explicit parallel(const dictionary &dict) : effectiveHeatConductivity(dict)
    {}

    // TypeName("parallel");

    tmp<volScalarField> kappaEff(std::span<kappaInfo> ki) const override
    {
        if (ki.size() <= 0) {
            return tmp<volScalarField>(nullptr);
        }
        auto alphaEff = max(ki[0].alpha, alpha_small_);
        tmp<volScalarField> ret = ki[0].kappa * alphaEff.ref();
        volScalarField alphaSum = alphaEff();
        volScalarField &result = ret.ref();
        for (size_t idx = 1; idx < ki.size(); idx++) {
            auto alphaEff = max(ki[idx].alpha, alpha_small_);
            result += ki[idx].kappa * alphaEff.ref();
            alphaSum += alphaEff.ref();
        }
        // alphaSum=max(alphaSum, this->alpha_small_);
        result /= alphaSum;
        return ret;
    }
};

class series : public effectiveHeatConductivity
{
  public:
    explicit series(const dictionary &dict) : effectiveHeatConductivity(dict) {}

    // TypeName("series");

    tmp<volScalarField> kappaEff(std::span<kappaInfo> ki) const override
    {
        if (ki.size() <= 0) {
            return tmp<volScalarField>(nullptr);
        }

        tmp<volScalarField> alpha_by_kappa_sum =
            max(ki[0].alpha, this->alpha_small_) / ki[0].kappa;
        volScalarField alphaSum = ki[0].alpha;
        for (size_t idx = 1; idx < ki.size(); idx++) {
            alpha_by_kappa_sum.ref() +=
                max(ki[idx].alpha, this->alpha_small_) / ki[idx].kappa;
            alphaSum += ki[idx].alpha;
        }
        alpha_by_kappa_sum.ref() = alphaSum / alpha_by_kappa_sum.ref();
        return alpha_by_kappa_sum;
    }
};

autoPtr<effectiveHeatConductivity> effectiveHeatConductivity::New(
    const dictionary &dict)
{
    using ret_t = autoPtr<effectiveHeatConductivity>;
    const word type = dict.lookup<word>("type");
    if (type == "parallel") {
        return ret_t(new parallel(dict));
    }
    if (type == "series") {
        return ret_t(new series(dict));
    }

    Info << "Unknown type " << type << " for effectiveHeatConducitivity"
         << endl;
    std::abort();

    return autoPtr<effectiveHeatConductivity>{nullptr};
}
