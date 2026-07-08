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

#include "capillary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

namespace Foam::functionObjects
{
defineTypeNameAndDebug(capillary, 0);
addToRunTimeSelectionTable(functionObject, capillary, dictionary);
} // namespace Foam::functionObjects

using namespace Foam;
using namespace Foam::functionObjects;

capillary::capillary(const word &name, const Time &runTime,
                     const dictionary &dict)
    : fvMeshFunctionObject{name, runTime, dict},
      write_capillary_pressure_{false}, write_phase_pressure_{false},
      write_capillary_force_{false}
{
    read(dict);
    Info << "Selecting capillary function" << endl;
}

bool capillary::read(const dictionary &dict)
{
    fvMeshFunctionObject::read(dict);
    write_capillary_pressure_ = dict.lookup<bool>("writeCapillaryPressure");
    write_phase_pressure_ = dict.lookup<bool>("writePhasePressure");
    write_capillary_force_ = dict.lookup<bool>("writeCapillaryForce");

    return true;
}

bool capillary::write()
{

    auto &cs =
        this->mesh().lookupObject<capillarySystem>(capillarySystem::dictName);
    auto &ps =
        this->mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName);
    auto &p = this->mesh().lookupObject<volScalarField>("p");

    for (auto &[to_phase, model] : cs.capillary_pressure_models()) {
        const volScalarField Pc("Pc_" +
                                    capillarySystem::phase_pair_name(
                                        cs.reference_phase(), to_phase),
                                model->capillary_pressure());
        if (write_phase_pressure_) {
            Pc.write();
        }
    }

    for (const auto &phase_name : cs.fluid_phase_names()) {

        if (write_phase_pressure_) {
            const volScalarField p_sub_pi =
                cs.pressure_difference_from_common_p(phase_name);
            const volScalarField pi("p_" + phase_name, p - p_sub_pi);
            pi.write();
        }

        if (write_capillary_force_) {
            const volVectorField force_i = cs.capillary_force(phase_name);
            force_i.write();
        }
    }
    return true;
}