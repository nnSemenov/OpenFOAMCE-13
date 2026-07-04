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

tmp<volScalarField> none::capillary_pressure() const
{
    auto &mesh = interface_.fluid().mesh();
    return volScalarField(
        IOobject("Pc" + interface_.name(), mesh.time().timePath(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar{dimPressure, scalar{0}});
}