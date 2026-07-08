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
    Foam::capillaryPressureModels::LeverettJ

Description
    Common capillary pressure model in porous media based on Leverett-J
    function. J is a function of liquid saturation (alpha_liq / porosity).
    Commonly used for gas, water/oil in rocks. Also, suitable of multiphase
    packed bed.

    Capillary pressure is calculated as follows:
        P_c = (1- porosity)/(porosity*dp) * sqrt(E1) * cos(theta)*J
        S_L = alpha_L / porosity
        J = J0 + beta * log((1-S_L)/S_L)

    This is not perfect. Currently specialized for packed bed. Will improve
    later.

    Supports multiple particle phases.

    Contact angle (contactAngle) is optional. If not assigned, cos(theta) is
    ignored. J is sometimes defined without this parameter.

    Example:
    \verbatim
    air_to_water {
        type LeverettJ;

        solids (solid);
        liquid water;

        J0 0.48;
        beta 0.036;

        E1 180;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: LeverettJ    | yes         |
        liquid       | Name of liquid phase    | yes         |
        solids       | List of solid phases    | yes         |
        contactAngle | contact angle in degree | no          |
        E1           | 1st ergun coefficient   | yes         |
        J0           | J function param        | yes         |
        beta         | J function param        | yes         |
    \endtable


SourceFiles
    LeverettJ.C

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_multiphaseEuler_LeverettJ_H
#define MIKENO_multiphaseEuler_LeverettJ_H

#include "capillaryPressureModel.H"

#include <optional>

namespace Foam::capillaryPressureModels
{
// Currently supports logit function only
class LeverettJ : public capillaryPressureModel
{
  protected:
    word liquidName_;
    List<word> solidNames_;

    scalar E1_;
    scalar J0_;
    scalar beta_;
    std::optional<scalar> contactAngle_;

  public:
    LeverettJ(const dictionary &dict, const phaseInterface &interface);
    virtual ~LeverettJ() override = default;

    TypeName("LeverettJ");

    [[nodiscard]] std::pair<tmp<volScalarField>, tmp<volScalarField>>
    capillary_pressure_with_derivative() const override;
};

} // namespace Foam::capillaryPressureModels

#endif // MIKENO_multiphaseEuler_LeverettJ_H
