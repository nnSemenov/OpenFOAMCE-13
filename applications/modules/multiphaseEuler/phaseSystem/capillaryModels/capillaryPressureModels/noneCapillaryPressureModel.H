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
    Foam::capillaryPressureModels::none

Description
    Fallback none capillary pressure model, always gives 0 pressure difference.

SourceFiles
    noneCapillaryPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_multiphaseEuler_nonCapillaryPressureModel_H
#define MIKENO_multiphaseEuler_nonCapillaryPressureModel_H

#include "capillaryPressureModel.H"

namespace Foam::capillaryPressureModels
{

class none : public capillaryPressureModel
{
  protected:
  public:
    TypeName("none");
    none(const dictionary &dict, const phaseInterface &interface);

    ~none() override = default;

    [[nodiscard]] tmp<volScalarField> capillary_pressure() const override;
};

} // namespace Foam::capillaryPressureModels

#endif // MIKENO_multiphaseEuler_nonCapillaryPressureModel_H
