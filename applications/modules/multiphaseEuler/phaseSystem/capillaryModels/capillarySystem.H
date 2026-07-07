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

Class
    Foam::capillarySystem

Description
    Capillary system organise pressure models in star-like. A phase is assigned
    to be reference phase. (only in this system, no relation to alpha solving),
    and all other fluid phases need to have a capillary model with it. The term
    `phaseA_to_phaseB` refers to capillary pressure $p_A - p_B$.

    Like most Euler-Euler multiphase solvers, `multiphaseEuler` assumes all
    phases share common pressure $p$. This is not breakable otherwise we rebuild
    everything. Here we use a convention to link common pressure with each
    fluid phase's own pressure: p = pressure of reference phase.

SourceFiles
    capillarySystem.C

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_multiphaseEuler_capillarySystem_H
#define MIKENO_multiphaseEuler_capillarySystem_H

#include <dictionary.H>
#include <phaseSystem.H>
#include <surfaceTensionCoefficientModel.H>

#include <memory>
#include <map>
#include <vector>

#include "capillaryPressureModel.H"

namespace Foam
{

class capillarySystem : public IOdictionary
{
  protected:
    const phaseSystem &phase_system_;

    // Reference phase for capillary.
    const word reference_phase_;
    std::map<word, autoPtr<capillaryPressureModels::capillaryPressureModel>>
        capillaryPressureModels_;

    void read_capillary_models(const dictionary &dict) &;

    void warn_if_movable_solids(const word &phase) const;

  public:
    explicit capillarySystem(const phaseSystem &phases);
    ~capillarySystem() override = default;

    TypeName("capillarySystem");

    inline const static word dictName{"capillaryProperties"};

    [[nodiscard]] inline auto &reference_phase() const &
    {
        return reference_phase_;
    }

    [[nodiscard]] inline auto &capillary_pressure_models() const & noexcept
    {
        return capillaryPressureModels_;
    }

    [[nodiscard]] inline static word phase_pair_name(const word &phaseA,
                                                     const word &phaseB)
    {
        return phaseA + "_to_" + phaseB;
    }

    [[nodiscard]] std::vector<word> fluid_phase_names() const;

    [[nodiscard]] std::vector<const phaseModel *> fluid_phases() const;

    // p - p_i. Difference between common pressure to phase pressure
    [[nodiscard]] volScalarField pressure_difference_from_common_p(
        const word &phase) const;

    [[nodiscard]] volScalarField pressure_difference_between_phases(
        const word &phase1, const word &phase2) const;

    [[nodiscard]] tmp<surfaceScalarField> Ff(const word &phase) const;

    [[nodiscard]] tmp<volVectorField> F(const word &phase) const;

    void add_to_Fs(PtrList<surfaceScalarField> &Fs) const;
    //- As Fs, but for the face-based algorithm
    void add_to_Ffs(PtrList<surfaceScalarField> &Ffs) const;

    [[nodiscard]] static autoPtr<capillarySystem> try_read(
        const phaseSystem &phaseSystem, const fvMesh &mesh);
};


} // namespace Foam

#endif
