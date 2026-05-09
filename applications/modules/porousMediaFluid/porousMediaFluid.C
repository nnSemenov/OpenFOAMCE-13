
#include <cstdlib>

#include "volFieldsFwd.H"
#include "porousMediaFluid.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace solvers {

    defineTypeNameAndDebug(porousMediaFluid, 0);
    addToRunTimeSelectionTable(solver, porousMediaFluid, fvMesh);
}
}

Foam::solvers::porousMediaFluid::porousMediaFluid(fvMesh&mesh) :
    fluid(mesh),
    U_physical(
        IOobject(
            "U_physical", 
            runTime.name(), 
            mesh,
            IOobject::READ_IF_PRESENT, 
            IOobject::AUTO_WRITE
        ),
        this->U_
    ),
    porosity_(
        IOobject(
            "porosity",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT, 
            IOobject::AUTO_WRITE
        ),
        mesh,
        scalar{1}
    ),
    porousPhases(mesh, runTime)
{
    const word fluid_phase_thermo_name = this->thermo().he().name();
    // Solid and fluid should have same energy type
    for(auto & porPhase: this->porousPhases.porousPhases()) {
        porPhase.second.thermo->validate("solid", fluid_phase_thermo_name);
    }

    this->updatePorosity();
    this->U_physical=this->U();
    this->U_physical /= this->porosity_;
}

Foam::solvers::porousMediaFluid::~porousMediaFluid() {}

void Foam::solvers::porousMediaFluid::updatePorosity() {
    // fill with 1
    porosity_.primitiveFieldRef()=1;
    porosity_.boundaryFieldRef()=1;
    for(auto & porPhase: porousPhases.porousPhases()) {
        porosity_-=porPhase.second.alpha;
    }
}


void Foam::solvers::porousMediaFluid::momentumPredictor() {
    isothermalFluid::momentumPredictor();
    this->U_physical=this->U();
    this->U_physical /= this->porosity_;
}
