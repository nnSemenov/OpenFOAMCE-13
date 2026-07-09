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
    Foam::fv::capillaryForce

Description
    Capillary force between vapor and liquids in porous media(rocks, packed bed,
    etc). Supports arbitrary fluid phases (but usually shouldn't be more than 4
    liquids, and physically very impossible for multiple gas phases).

Example
    capillaryForce {
        type capillaryForce;
        libs ("libmultiphaseEulerFvModels.so");

        adjustTimeStep yes;
        maxFo 1;
    }

SourceFiles
    capillaryForce.C

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_multiphaseEuler_capillaryForce_H
#define MIKENO_multiphaseEuler_capillaryForce_H

#include <fvModel.H>
#include <optional>

#include "capillarySystem.H"

namespace Foam::fv
{
class capillaryForce : public fvSource
{
  protected:
    std::optional<scalar> maxFo_;

    void read_coeffs(const dictionary &dict) &;

  public:
    capillaryForce(const word &name, const word &modelType, const fvMesh &mesh,
                   const dictionary &dict);
    ~capillaryForce() override = default;
    TypeName("capillaryForce");

    [[nodiscard]] const capillarySystem &get_capillary_system() const;

    [[nodiscard]] word match_phase_from_field(const capillarySystem &,
                                              const word &fieldName) const;

    [[nodiscard]] bool addsSupToField(const word &fieldName) const override;

    void addSup(const volScalarField &alphai, const volScalarField &rhoi,
                const volVectorField &Ui, fvVectorMatrix &UiEqn) const override;

    inline bool movePoints() override { return false; };
    inline void topoChange(const polyTopoChangeMap &) override {};
    inline void mapMesh(const polyMeshMap &) override {};
    inline void distribute(const polyDistributionMap &) override {};

    bool read(const dictionary &dict) override;

    [[nodiscard]] scalar maxDeltaT() const override;
};
} // namespace Foam::fv

#endif // MIKENO_multiphaseEuler_capillaryForce_H
