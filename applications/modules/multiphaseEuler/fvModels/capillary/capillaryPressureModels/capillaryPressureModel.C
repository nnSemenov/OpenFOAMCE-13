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


#include "capillaryPressureModel.H"
#include <phaseInterface.H>

namespace Foam::capillaryPressureModels
{
defineTypeNameAndDebug(capillaryPressureModel, 0);
defineRunTimeSelectionTable(capillaryPressureModel, dictionary);
} // namespace Foam::capillaryPressureModels

Foam::autoPtr<Foam::capillaryPressureModels::capillaryPressureModel> Foam::
    capillaryPressureModels::capillaryPressureModel::New(
        const dictionary &dict, const phaseInterface &interface)
{
    const word type = dict.lookup("type");
    Info << indentOrNl << "Selecting capillaryPressureModel for "
         << interface.name() << ": " << type << endl;

    auto iter = dictionaryConstructorTablePtr_->find(type);
    if (iter == dictionaryConstructorTablePtr_->end()) {
        FatalIOErrorInFunction(dict)
            << "Unknown capillaryPressureModel type " << type << endl
            << endl
            << "Valid capillaryPressureModel types are : "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return iter()(dict, interface);
}

Foam::capillaryPressureModels::capillaryPressureModel::capillaryPressureModel(
    const dictionary &dict, const phaseInterface &interface)
    : interface_{interface}, sigma_model_{sigmaModel(dict, interface)}
{}

const Foam::surfaceTensionCoefficientModel &Foam::capillaryPressureModels::
    capillaryPressureModel::sigmaModel(const dictionary &dict,
                                       const phaseInterface &interface)
{

    auto &table = interface.fluid().surfaceTensionCoefficientModels();
    const phaseInterfaceKey key{interface.phase1(), interface.phase2()};

    if (not table.found(key)) {
        FatalIOErrorInFunction(dict)
            << "Missing surface tension between phase " << key.first()
            << " and " << key.second() << exit(FatalIOError);
    }

    return table[key];
}