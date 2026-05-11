#include <cstdio>
#include <cstddef>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <rapidcsv.h>

#include "argList.H"

#include "fluidThermo.H"
#include "fluidMulticomponentThermo.H"
#include "compressibleMomentumTransportModel.H"

#include "specie.H"

using Foam::List;
using Foam::scalar;
using Foam::scalarList;

struct testData
{
    std::vector<float> p;
    std::vector<float> T;
    std::vector<float> rho;
    std::vector<float> Cpv;
    // Optional
    std::vector<std::vector<float>> Y_;
    std::vector<float> mu_;
    std::vector<float> kappa_;

    ~testData() = default;

    void delete_data(size_t index) &
    {
        this->p.erase(this->p.begin() + index);
        this->T.erase(this->T.begin() + index);
        this->rho.erase(this->rho.begin() + index);
        this->Cpv.erase(this->Cpv.begin() + index);

        for (auto &vec : Y_) {
            vec.erase(vec.begin() + index);
        }
        if (with_mu()) {
            mu_.erase(this->mu_.begin() + index);
        }
        if (with_kappa()) {
            kappa_.erase(this->kappa_.begin() + index);
        }
    }

    [[nodiscard]] bool mutlicomponent() const { return not Y_.empty(); }

    [[nodiscard]] bool with_mu() const { return not mu_.empty(); }

    [[nodiscard]] bool with_kappa() const { return not kappa_.empty(); }
};

struct dataset_load_option
{
    const Foam::speciesTable *species_table{nullptr};
    std::string specific_heat_name{"CP"};
    bool read_transport{false};
};

std::unique_ptr<testData> load_data(const std::string &filename,
                                    const dataset_load_option &option);

scalar relativeDiff(const scalar value, const scalar accurate)
{
    return (value - accurate) / accurate;
}

struct relativeDiffRange
{
    scalar min{0};
    scalar max{0};

    void update(scalar rDiff)
    {
        this->min = Foam::min(this->min, rDiff);
        this->max = Foam::max(this->max, rDiff);
    }
};

int main(int argc, char *argv[])
{
    using namespace Foam;

    argList::addOption("eos", "name", "Equation of state for thermodynamic");
    argList::addBoolOption("mixture", "If is multicomponent");
    argList::addBoolOption("transport", "Read and validate mu and kappa");

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    autoPtr<fluidThermo> thermoPtr_;
    const speciesTable *specie_table = nullptr;
    if (args.optionFound("mixture")) {
        thermoPtr_ = fluidMulticomponentThermo::New(mesh).ptr();
        const auto &sp_table =
            dynamic_cast<fluidMulticomponentThermo &>(thermoPtr_()).species();
        std::printf("Specie list: [");
        forAll(sp_table, spIdx) { std::printf("%s ", sp_table[spIdx].c_str()); }
        std::printf("\b]\n");
        specie_table =
            &dynamic_cast<fluidMulticomponentThermo &>(thermoPtr_()).species();
    }
    else {
        thermoPtr_ = fluidThermo::New(mesh);
    }

    const word eos = args.option("eos");
    dataset_load_option load_opt{
        .species_table = specie_table,
        .read_transport = args.optionFound("transport"),
    };
    if (thermoPtr_->he().name() == "e") {
        load_opt.specific_heat_name = "CV";
    }
    const auto data = load_data("data/" + eos + ".csv", load_opt);
    if (data == nullptr) {
        return 1;
    }
    if ((specie_table not_eq nullptr) not_eq data->mutlicomponent()) {
        std::fprintf(stderr,
                     "Fatal error: data and thermodynamics doesn't match\n "
                     "Must be either multicomponent or pure.");
        return 1;
    }

    auto &p = thermoPtr_->p();
    auto &T = thermoPtr_->T();
    auto &Cpv = thermoPtr_->Cpv();
    auto &mu = thermoPtr_->mu();
    auto &kappa = thermoPtr_->kappa();

    if (data->T.size() > mesh.nCells()) {
        Info << "No enough cells: required " << data->T.size()
             << ", but only have " << mesh.nCells() << endl;
        return 1;
    }
    forAll(data->T, idx)
    {
        p[idx] = data->p[idx];
        T[idx] = data->T[idx];
    }
    if (data->mutlicomponent()) {
        assert(specie_table);
        const auto &Y = data.get()->Y_;
        auto &thermo = dynamic_cast<fluidMulticomponentThermo &>(thermoPtr_());
        auto &Y_dest = thermo.Y();
        const auto &spTable = *specie_table;
        // const auto & spTable = thermo.species();
        forAll(spTable, sp_index)
        {
            forAll(data->T, idx) { Y_dest[sp_index][idx] = Y[sp_index][idx]; }
        }
        thermo.normaliseY();
    }
    // correct thermo
    thermoPtr_->he() = thermoPtr_->he(p, T);
    thermoPtr_->correct();

    const auto rho = thermoPtr_->rho();

    relativeDiffRange rhoDiff, CpvDiff;
    relativeDiffRange muDiff, kappaDiff;

    forAll(data->T, cell)
    {
        std::printf("p = %.1e Pa, T = %.2f K, ", p[cell], T[cell]);

        if (data->mutlicomponent()) {
            auto &thermo =
                dynamic_cast<fluidMulticomponentThermo &>(thermoPtr_());
            auto &Y_dest = thermo.Y();
            forAll(*specie_table, spidx)
            {
                std::printf("Y_%s = %1.2e, ", (*specie_table)[spidx].c_str(),
                            Y_dest[spidx][cell]);
            }
        }

        const scalar rhoRDiff = relativeDiff(rho()[cell], data->rho[cell]);
        rhoDiff.update(rhoRDiff);
        std::printf("rho diff = %1.2e, ", rhoRDiff);

        const scalar CpvRelDiff = relativeDiff(Cpv[cell], data->Cpv[cell]);
        CpvDiff.update(CpvRelDiff);
        std::printf("Cpv diff = %1.2e, ", CpvRelDiff);

        if (load_opt.read_transport) {
            const scalar muRDiff = relativeDiff(mu[cell], data->mu_[cell]);
            muDiff.update(muRDiff);
            std::printf("mu diff = %1.2e, ", muRDiff);
            const scalar kappaRDiff =
                relativeDiff(kappa[cell], data->kappa_[cell]);
            kappaDiff.update(kappaRDiff);
            std::printf("kappa diff = %1.2e, ", kappaRDiff);
        }

        std::printf("\n");
    }
    std::printf("\n\n\nSummary:\n");
    std::printf("\tRelative diff of rho: [%1.2e, %1.2e]\n", rhoDiff.min,
                rhoDiff.max);
    std::printf("\tRelative diff of Cpv: [%1.2e, %1.2e]\n", CpvDiff.min,
                CpvDiff.max);
    if (load_opt.read_transport) {
        std::printf("\tRelative diff of mu: [%1.2e, %1.2e]\n", muDiff.min,
                    muDiff.max);
        std::printf("\tRelative diff of kappa: [%1.2e, %1.2e]\n", kappaDiff.min,
                    kappaDiff.max);
    }
    //    Info<<endl;
    return 0;
}

std::unique_ptr<testData> load_data(const std::string &filename,
                                    const dataset_load_option &option)
{
    std::unique_ptr<testData> ret{nullptr};
    const auto species_table = option.species_table;
    const bool contains_specie = species_table;
    try {
        rapidcsv::Document csv(filename);
        ret = std::make_unique<testData>();

        ret->T = csv.GetColumn<float>("T/K");
        ret->p = csv.GetColumn<float>("p/Pa");
        ret->Cpv = csv.GetColumn<float>(option.specific_heat_name + "/J/kg/K");
        ret->rho = csv.GetColumn<float>("rho/kg/cum");

        const size_t N_data = ret->T.size();

        if (contains_specie) {
            auto &data = *ret;
            const size_t N_specie = species_table->size();
            data.Y_.resize(species_table->size());
            forAll(*species_table, specie_index)
            {
                std::string col_name = "Y_";
                col_name += (*species_table)[specie_index].c_str();
                data.Y_[specie_index] = csv.GetColumn<float>(col_name);
            }
            // normalize data
            for (size_t idx = 0; idx < N_data; idx++) {
                scalar Ysum = 0;
                for (size_t spIdx = 0; spIdx < N_specie; spIdx++) {
                    Ysum += data.Y_[spIdx][idx];
                }
                if (not(Ysum > 0)) {
                    throw std::runtime_error("Found non-positive Y sum = " +
                                             std::to_string(Ysum));
                }

                for (size_t spIdx = 0; spIdx < N_specie; spIdx++) {
                    data.Y_[spIdx][idx] /= Ysum;
                }
            }
        }

        if (option.read_transport) {
            ret->mu_ = csv.GetColumn<float>("MU/Pa-sec");
            ret->kappa_ = csv.GetColumn<float>("K/W/m/K");
        }

        std::vector<size_t> erase_idx;
        for (size_t idx = 0; idx < N_data; idx++) {
            if (not(ret->rho[idx] > 0) or not(ret->Cpv[idx] > 0)) {
                erase_idx.emplace_back(idx);
                continue;
            }
            if (ret->T[idx] > 5000) {
                erase_idx.emplace_back(idx);
                continue;
            }
        }

        for (auto it = erase_idx.rbegin(); it not_eq erase_idx.rend(); ++it) {
            ret->delete_data(*it);
        }

        std::printf("Erased %zu from %zu lines\n", erase_idx.size(), N_data);
    }
    catch (const std::exception &e) {
        std::fprintf(stderr, "Exception while readinig %s: %s\n",
                     filename.c_str(), e.what());
        return nullptr;
    }

    return ret;
}
