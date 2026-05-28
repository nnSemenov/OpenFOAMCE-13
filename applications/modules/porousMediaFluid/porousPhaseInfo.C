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

#include "porousPhaseInfo.H"
#include <algorithm>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <variant>
#include <vector>

#include "cubicEOS.H"

using Foam::word, Foam::List;

const Foam::word Foam::porousPhaseInfo::FLUID_PHASE_NAME{"fluid"};

Foam::heatTransfer Foam::parse_heatTransfer(const Foam::dictionary &dict)
{
    using namespace Foam;
    const word type = dict.lookup<word>("type");
    if (type == "equilibrium") {
        return equilibriumHeatTransfer{};
    }

    if (type == "nonEquilibrium") {
        const auto alpha = dict.lookup<dimensionedScalar>("heatTransferCoeff");
        const auto area = dict.lookup<dimensionedScalar>("specificSurfaceArea");
        // Info<<"Dimension of heat transfer coeff: "<<alpha.dimensions()<<"\n,
        // specific surface area: "<<area.dimensions()<<endl;
        const dimensionSet dimHeatTransferCoeff(1, 0, -3, -1, 0);
        const dimensionSet dimSpecificSurfaceArea(0, -1, 0, 0, 0);
        bool ok = true;
        if (not(alpha.dimensions() == dimHeatTransferCoeff)) {
            Info << "Wrong dimension for heatTransferCoeff, must be "
                 << dimHeatTransferCoeff << ", found " << alpha.dimensions()
                 << endl;
            ok = false;
        }
        if (not(area.dimensions() == dimSpecificSurfaceArea)) {
            Info << "Wrong dimension for specificSurfaceArea, must be "
                 << dimSpecificSurfaceArea << ", found " << area.dimensions()
                 << endl;
            ok = false;
        }
        if (not ok) {
            std::abort();
        }

        const bool explicitTerm = dict.lookupOrDefault<bool>("explicit", false);

        return nonEquilibriumHeatTransfer{.heatTransferCoefficient = alpha,
                                          .effectiveSpecificSurfaceArea = area,
                                          .explicitTerm = explicitTerm};
    }

    Info << "Unknown heat transfer type " << type << endl;
    std::abort();
    return {};
}

Foam::porousPhaseInfo::porousPhaseInfo(fvMesh &mesh, const Time &runTime)
    : heatTransferGraph_(0)
{
    IOdictionary porousPhaseDict(IOobject("phaseProperties", runTime.constant(),
                                          mesh, IOobject::MUST_READ,
                                          IOobject::NO_WRITE));
    wordList phaseNameList(porousPhaseDict.lookup<wordList>("porousPhases"));

    for (const word &phaseName : phaseNameList) {
        if (phaseName == FLUID_PHASE_NAME) {
            FatalErrorInFunction
                << "Invalid phase name for porous media: " << phaseName << endl;
            ::Foam::abort(::Foam::FatalError);
        }

        if (porousPhases_.find(phaseName) not_eq porousPhases_.end()) {
            FatalErrorInFunction << "Porous phase name " << phaseName
                                 << " is already defined." << endl;
            ::Foam::abort(::Foam::FatalError);
        }
        auto thermo = solidThermo::New(mesh, phaseName);
        auto transport = solidThermophysicalTransportModel::New(thermo());
        thermo->validate("solid", "h", "e");
        auto it = porousPhases_.emplace(
            phaseName,
            porousPhase{.alpha = volScalarField(
                            IOobject("alpha." + phaseName, runTime.name(), mesh,
                                     IOobject::MUST_READ, IOobject::AUTO_WRITE,
                                     true // registered
                                     ),
                            mesh),
                        .thermo = thermo,
                        .thermophysicalTransport = transport});
        // force solver to write alpha (why AUTO_WRITE above doesn't work ?!)
        it.first->second.alpha.writeOpt() = IOobject::AUTO_WRITE;
        const word pN = it.first->second.thermo->phaseName();
        assert(pN == phaseName);
    }

    Info << "Create " << phaseNameList.size()
         << " porous phases: " << phaseNameList << endl;

    for (auto &phase : this->porousPhases_) {
        auto &g = this->heatTransferGraph_;
        auto vert = boost::add_vertex(phase.first, g);
        phase.second.vertex = vert;
    }
    this->fluid_vertex =
        boost::add_vertex(FLUID_PHASE_NAME, this->heatTransferGraph_);

    auto getPhaseVertex = [this](const word &phaseName) {
        if (phaseName == FLUID_PHASE_NAME) {
            return this->fluid_vertex;
        }
        return this->porousPhases_.at(phaseName).vertex;
    };
    {
        const auto &heatTransferDict = porousPhaseDict.subDict("heatTransfer");
        const auto keys = heatTransferDict.keys(true);
        for (const word &key : keys) {
            word phaseA, phaseB;
            const word error = parseBinarySpeciePair(key, phaseA, phaseB);
            if (not error.empty()) {
                FatalErrorInFunction
                    << "Invalid phase heat transfer table: " << error << endl;
                ::Foam::abort(::Foam::FatalError);
            }
            Info << "Creating heat transfer between " << phaseA << " and "
                 << phaseB << endl;
            for (const auto &phaseName : {phaseA, phaseB}) {
                if (not isValidPhase(phaseName)) {
                    FatalErrorInFunction << "Invalid phase name \"" << phaseName
                                         << "\", no such phase" << endl;
                    ::Foam::abort(::Foam::FatalError);
                }
            }

            auto heatTransfer =
                parse_heatTransfer(heatTransferDict.subDict(key));

            auto &g = this->heatTransferGraph_;
            boost::add_edge(
                getPhaseVertex(phaseA), getPhaseVertex(phaseB),
                phaseTransferProperty{.heatTransfer_ = heatTransfer}, g);
        }
    }
    Info << "Added " << boost::num_edges(this->heatTransferGraph_)
         << " heat transfer pathes" << endl;

    // Create effective kappa models
    {
        const dictionary &eKDict =
            porousPhaseDict.subDict("effectiveHeatConductivity");
        List<keyType> keys = eKDict.keys(true);
        keys.append(eKDict.keys(false));
        for (const word &key : keys) {
            const dictionary &dict = eKDict.subDict(key);
            auto ek = effectiveHeatConductivity::New(dict);
            this->effectiveHeatConductivityModels.emplace(key, ek);
        }
        Info << "Loaded " << keys.size()
             << " effective heat conductivity model(s): \n"
             << keys << endl;
    }

    auto summary = this->heatTransferSummary();
    Info << "Heat transfer summary: \n";
    for (auto &node : summary) {
        Info << "  (";
        for (const word &phaseName : node.thermalEquilibriumPhaseNames) {
            Info << phaseName << ", ";
        }
        Info << "\b\b) has following interaction(s): \n";
        for (const auto &source : node.heatTransferSources) {
            Info << "    "
                 << node.thermalEquilibriumPhaseNames[source.phaseA_index]
                 << "<->" << source.phaseB_name << '\n';
        }
        Info << '\n';
    }
    Info << endl;
}

class Visitor : public boost::default_bfs_visitor
{
  public:
    std::vector<word> *phase_list;
    std::vector<
        Foam::porousPhaseInfo::heatTransferGraph_type::vertex_descriptor>
        *vert_list;

    template <typename vertex, typename Graph>
    void discover_vertex(vertex v, Graph g) const
    {
        word phaseName = g[v];
        phase_list->emplace_back(phaseName);
        vert_list->emplace_back(v);
    }
};

std::vector<std::vector<word>> Foam::porousPhaseInfo::thermalEquilibriumPhases()
    const
{
    heatTransferGraph_type graph{this->heatTransferGraph_};

    // Remove all non-thermal-equilibrium heat transfer links
    auto pred = [&](auto edge) -> bool {
        const heatTransfer &ht = graph[edge].heatTransfer_;
        if (std::get_if<equilibriumHeatTransfer>(&ht) == nullptr) {
            return true;
        }
        return false;
    };
    boost::remove_edge_if(pred, graph);

    std::vector<std::vector<word>> ret;
    while (true) {
        if (boost::num_vertices(graph) <= 0) {
            break;
        }
        std::vector<word> phases;
        std::vector<heatTransferGraph_type::vertex_descriptor> verts;
        Visitor vis;
        vis.phase_list = &phases;
        vis.vert_list = &verts;
        heatTransferGraph_type::vertex_iterator start_node =
            boost::vertices(graph).first;
        boost::breadth_first_search(graph, *start_node, boost::visitor(vis));

        std::sort(verts.begin(), verts.end());
        for (auto it = verts.rbegin(); it not_eq verts.rend(); ++it) {
            boost::remove_vertex(*it, graph);
        }
        assert(phases.size() > 0);
        std::sort(phases.begin(), phases.end());
        ret.emplace_back(std::move(phases));
    }
    assert(ret.size() > 0);

    return ret;
}

Foam::porousPhaseInfo::heatTransferGraph_type::vertex_descriptor Foam::
    porousPhaseInfo::vertexOfPhase(const word &phaseName) const
{
    if (phaseName == FLUID_PHASE_NAME) {
        return fluid_vertex;
    }
    return porousPhases_.at(phaseName).vertex;
}

const Foam::heatTransfer *Foam::porousPhaseInfo::heatTransferBetween(
    const word &phaseA, const word &phaseB) const
{
    auto vertA = this->vertexOfPhase(phaseA);
    auto vertB = this->vertexOfPhase(phaseB);

    auto pair = boost::edge(vertA, vertB, heatTransferGraph_);
    if (pair.second) { // edge exist
        return &heatTransferGraph_[pair.first].heatTransfer_;
    }
    return nullptr;
}

const Foam::nonEquilibriumHeatTransfer *Foam::porousPhaseInfo::
    nonEquilibriumHeatTransferBetween(const word &phaseA,
                                      const word &phaseB) const
{
    const Foam::heatTransfer *ht = heatTransferBetween(phaseA, phaseB);
    return std::get_if<nonEquilibriumHeatTransfer>(ht);
}

Foam::word Foam::porousPhaseInfo::combinedPhaseName(
    std::vector<word> phaseNames)
{
    std::sort(phaseNames.begin(), phaseNames.end());
    word phaseList;
    for (const word &ph : phaseNames) {
        phaseList += ph;
        phaseList += ':';
    }
    phaseList.pop_back();
    return phaseList;
}

Foam::effectiveHeatConductivity *Foam::porousPhaseInfo::
    effectiveKappaModelForPhaseList(std::vector<word> phaseNames) const
{
    const word combinedName = combinedPhaseName(phaseNames);

    auto it = this->effectiveHeatConductivityModels.find(combinedName);
    if (it not_eq this->effectiveHeatConductivityModels.end()) {
        return const_cast<Foam::effectiveHeatConductivity *>(&it->second());
    }
    it = this->effectiveHeatConductivityModels.find("default");
    if (it not_eq this->effectiveHeatConductivityModels.end()) {
        return const_cast<Foam::effectiveHeatConductivity *>(&it->second());
    }
    return nullptr;
}

std::vector<Foam::porousPhaseInfo::heatTransferInfoNode> Foam::porousPhaseInfo::
    heatTransferSummary() const
{
    std::vector<Foam::porousPhaseInfo::heatTransferInfoNode> ret;
    {
        auto eq_phase_list = this->thermalEquilibriumPhases();
        ret.reserve(eq_phase_list.size());

        for (auto &phaseNames : eq_phase_list) {
            effectiveHeatConductivity *effKappaModel =
                this->effectiveKappaModelForPhaseList(phaseNames);
            if (phaseNames.size() > 1 and effKappaModel == nullptr) {
                const word key = this->combinedPhaseName(phaseNames);
                Info << "Missing effective heat conducitivty model for "
                        "multiple thermal equilibrium phases: "
                     << key;
                Info << "\n You should specify effective heat conducitivty "
                        "model either for this phase combination or \"default\""
                     << endl;
                std::abort();
            }
            ret.emplace_back(heatTransferInfoNode{
                .thermalEquilibriumPhaseNames = std::move(phaseNames),
                .effectiveKappaModel = effKappaModel,
                .heatTransferSources = {}});
        }
    }
    const size_t N = ret.size();
    for (size_t i = 0; i < N; i++) {
        for (size_t phaseA_index = 0;
             phaseA_index < ret[i].thermalEquilibriumPhaseNames.size();
             phaseA_index++) {
            const word &phaseA =
                ret[i].thermalEquilibriumPhaseNames[phaseA_index];
            for (size_t j = 0; j < i; j++) {
                for (size_t phaseB_index = 0;
                     phaseB_index < ret[j].thermalEquilibriumPhaseNames.size();
                     phaseB_index++) {
                    const word &phaseB =
                        ret[j].thermalEquilibriumPhaseNames[phaseB_index];
                    const nonEquilibriumHeatTransfer *nonEq =
                        nonEquilibriumHeatTransferBetween(phaseA, phaseB);
                    if (nonEq == nullptr) {
                        continue;
                    }
                    ret[i].heatTransferSources.emplace_back(phaseA_index,
                                                            phaseB, *nonEq);
                    ret[j].heatTransferSources.emplace_back(phaseB_index,
                                                            phaseA, *nonEq);
                }
            }
        }
    }

    return ret;
}