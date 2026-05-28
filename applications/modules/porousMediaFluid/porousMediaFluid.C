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

#include <cstdlib>

#include "volFieldsFwd.H"
#include "porousMediaFluid.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace solvers
{

defineTypeNameAndDebug(porousMediaFluid, 0);
addToRunTimeSelectionTable(solver, porousMediaFluid, fvMesh);
} // namespace solvers
} // namespace Foam

Foam::solvers::porousMediaFluid::porousMediaFluid(fvMesh &mesh)
    : fluid(mesh),
      U_physical(IOobject("U_physical", runTime.name(), mesh,
                          IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                 this->U_),
      porosity_(IOobject("porosity", runTime.name(), mesh,
                         IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                mesh, scalar{1}),
      porousPhases(mesh, runTime)
{
    const word fluid_phase_thermo_name = this->thermo().he().name();
    // Solid and fluid should have same energy type
    for (auto &porPhase : this->porousPhases.porousPhases()) {
        porPhase.second.thermo->validate("solid", fluid_phase_thermo_name);
    }

    this->updatePorosity();
    this->U_physical = this->U();
    this->U_physical /= this->porosity_;
}

Foam::solvers::porousMediaFluid::~porousMediaFluid() {}

void Foam::solvers::porousMediaFluid::updatePorosity()
{
    // fill with 1
    porosity_.primitiveFieldRef() = 1;
    porosity_.boundaryFieldRef() = 1;
    for (auto &porPhase : porousPhases.porousPhases()) {
        porosity_ -= porPhase.second.alpha;
    }
}

void Foam::solvers::porousMediaFluid::momentumPredictor()
{
    isothermalFluid::momentumPredictor();
    this->U_physical = this->U();
    this->U_physical /= this->porosity_;
}
