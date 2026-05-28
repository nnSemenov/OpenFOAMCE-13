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

#ifndef porousPhaseInfo_H
#define porousPhaseInfo_H

#include <variant>
#include <map>
#include <optional>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_graph.hpp>

#include "effectiveHeatConductivity.H"

#include "volFieldsFwd.H"
#include "IOdictionary.H"
#include "solidThermophysicalTransportModel.H"

namespace Foam
{

struct equilibriumHeatTransfer
{};
struct nonEquilibriumHeatTransfer
{
    dimensionedScalar heatTransferCoefficient;
    dimensionedScalar effectiveSpecificSurfaceArea;
    bool explicitTerm;
};
using heatTransfer =
    std::variant<equilibriumHeatTransfer, nonEquilibriumHeatTransfer>;

heatTransfer parse_heatTransfer(const dictionary &dict);

class porousPhaseInfo
{
  public:
    static const word FLUID_PHASE_NAME;

    struct phaseTransferProperty
    {
        heatTransfer heatTransfer_;
    };

    using heatTransferGraph_type =
        boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              std::string, phaseTransferProperty,
                              boost::no_property>;

    struct porousPhase
    {
        volScalarField alpha;
        autoPtr<solidThermo> thermo;
        autoPtr<solidThermophysicalTransportModel> thermophysicalTransport;
        heatTransferGraph_type::vertex_descriptor vertex{};
    };

  protected:
    std::map<word, porousPhase> porousPhases_;

    heatTransferGraph_type heatTransferGraph_;
    heatTransferGraph_type::vertex_descriptor fluid_vertex{};

    std::map<word, autoPtr<effectiveHeatConductivity>>
        effectiveHeatConductivityModels;

  public:
    explicit porousPhaseInfo(fvMesh &mesh, const Time &runTime);

    porousPhaseInfo(const porousPhaseInfo &) = delete;

    ~porousPhaseInfo() = default;

    auto &porousPhases() & { return porousPhases_; }

    const auto &porousPhases() const & { return porousPhases_; }

    bool isValidPhase(const word &phaseName) const
    {
        return this->porousPhases_.contains(phaseName) or
            (phaseName == FLUID_PHASE_NAME);
    }
    //- List of thermal equilibrium phases
    std::vector<std::vector<word>> thermalEquilibriumPhases() const;

    heatTransferGraph_type::vertex_descriptor vertexOfPhase(
        const word &phaseName) const;

    const heatTransfer *heatTransferBetween(const word &phaseA,
                                            const word &phaseB) const;

    const nonEquilibriumHeatTransfer *nonEquilibriumHeatTransferBetween(
        const word &phaseA, const word &phaseB) const;

    static word combinedPhaseName(std::vector<word> phaseNames);

    effectiveHeatConductivity *effectiveKappaModelForPhaseList(
        std::vector<word> phaseNames) const;

    struct heatTransferInfoNode
    {
        struct heatTransferSource
        {
            size_t phaseA_index;
            word phaseB_name;
            const nonEquilibriumHeatTransfer &heatTransferInfo;
        };

        std::vector<word> thermalEquilibriumPhaseNames;
        // Not used for single phase
        effectiveHeatConductivity *effectiveKappaModel;
        std::vector<heatTransferSource> heatTransferSources;

        std::optional<size_t> fluidPhaseIndex() const
        {
            for (size_t idx = 0; idx < thermalEquilibriumPhaseNames.size();
                 idx++) {
                const word &phaseName = thermalEquilibriumPhaseNames[idx];
                if (phaseName == FLUID_PHASE_NAME) {
                    return idx;
                }
            }
            return std::nullopt;
        }
    };
    std::vector<heatTransferInfoNode> heatTransferSummary() const;
};

} // namespace Foam

#endif // porousPhaseInfo_H