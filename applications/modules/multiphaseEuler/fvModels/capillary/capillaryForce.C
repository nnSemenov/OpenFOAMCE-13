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


#include "capillaryForce.H"
#include <addToRunTimeSelectionTable.H>
#include <fvcGrad.H>

namespace Foam::fv
{
defineTypeNameAndDebug(capillaryForce, 0);
addToRunTimeSelectionTable(fvModel, capillaryForce, dictionary);
} // namespace Foam::fv

using namespace Foam;

fv::capillaryForce::capillaryForce(const word &name, const word &modelType,
                                   const fvMesh &mesh, const dictionary &dict)
    : fvSource{name, modelType, mesh, dict}, capillarySystem_{nullptr}
{
    this->read_capillary_model(dict);
}

void fv::capillaryForce::read_capillary_model(const dictionary &dict) &
{
    this->capillarySystem_ = new capillarySystem{
        this->mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName),
        dict};
}

bool fv::capillaryForce::read(const dictionary &dict)
{
    if (fvSource::read(dict)) {
        read_capillary_model(dict);
        return true;
    }
    return false;
}

word fv::capillaryForce::match_phase_from_field(const word &fieldName) const
{

    auto phase_names = this->capillarySystem_->fluid_phase_names();
    for (auto &name : phase_names) {
        if (fieldName == IOobject::groupName("U", name)) {
            return name;
        }
    }
    return word::null;
}

bool fv::capillaryForce::addsSupToField(const word &fieldName) const
{
    const auto matched_phase = match_phase_from_field(fieldName);
    return not matched_phase.empty();
}

void fv::capillaryForce::addSup(const volScalarField &alphai,
                                const volScalarField &rhoi,
                                const volVectorField &Ui,
                                fvVectorMatrix &UiEqn) const
{
    const word phase_name = match_phase_from_field(Ui.name());

    if (phase_name.empty()) {
        FatalErrorInFunction
            << "Trying to apply capillary force on field " << Ui.name()
            << " but it's not a fluid velocity field" << exit(FatalError);
        return;
    }
    if (alphai.name() not_eq IOobject::groupName("alpha", phase_name)) {
        FatalErrorInFunction << "Vol-fraction field " << alphai.name()
                             << " shouldn't be with velocity field "
                             << Ui.name() << exit(FatalError);
        return;
    }
    if (rhoi.name() not_eq IOobject::groupName("rho", phase_name)) {
        FatalErrorInFunction << "Density field " << rhoi.name()
                             << " shouldn't be with velocity field "
                             << Ui.name() << exit(FatalError);
        return;
    }

    auto deltaP =
        capillarySystem_->pressure_difference_from_common_p(phase_name);

    auto force = alphai * fvc::grad(deltaP, "grad(Pc)");
    UiEqn += force;
}