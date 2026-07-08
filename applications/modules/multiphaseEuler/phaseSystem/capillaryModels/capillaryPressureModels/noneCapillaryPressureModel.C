/*---------------------------------------------------------------------------*\
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


#include "noneCapillaryPressureModel.H"

#include <addToRunTimeSelectionTable.H>

namespace Foam::capillaryPressureModels
{
defineTypeNameAndDebug(none, 0);
addToRunTimeSelectionTable(capillaryPressureModel, none, dictionary);
} // namespace Foam::capillaryPressureModels

using namespace Foam;
using namespace capillaryPressureModels;

none::none(const dictionary &dict, const phaseInterface &interface)
    : capillaryPressureModel{dict, interface}
{}

std::pair<tmp<volScalarField>, tmp<volScalarField>> none::
    capillary_pressure_with_derivative() const
{
    auto &mesh = interface_.fluid().mesh();

    return std::make_pair(
        volScalarField::New("Pc" + interface_.name(), mesh,
                            dimensionedScalar{dimPressure, Zero}),
        volScalarField ::New(Pc_derivative_name + interface_.name(), mesh,
                             dimensionedScalar{dimPressure, Zero}));
}