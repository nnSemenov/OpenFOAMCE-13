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
    Foam::capillaryPressureModels::capillaryPressureModel

Description
    Base class of physical models for capillary pressure between 2 fluid phases.

SourceFiles
    capillaryPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_multiphaseEuler_capillaryPressureModel_H
#define MIKENO_multiphaseEuler_capillaryPressureModel_H

#include <dictionary.H>
#include <phaseSystem.H>
#include <surfaceTensionCoefficientModel.H>
#include <runTimeSelectionTables.H>

#include <utility>

namespace Foam::capillaryPressureModels
{

class capillaryPressureModel
{
  protected:
    const phaseInterface interface_;
    const surfaceTensionCoefficientModel &sigma_model_;

  public:
    TypeName("capillaryPressureSystem");
    //    capillaryPressureModel(const phaseModel &phaseA, const phaseModel
    //    &phaseB,
    //                           const phaseSystem &ps, const dictionary &dict);

    capillaryPressureModel(const dictionary &dict,
                           const phaseInterface &interface);

    virtual ~capillaryPressureModel() = default;

    static const surfaceTensionCoefficientModel &sigmaModel(
        const dictionary &dict, const phaseInterface &interface);

    [[nodiscard]] inline auto sigma() const { return sigma_model_.sigma(); }

    [[nodiscard]] inline tmp<volScalarField> capillary_pressure() const
    {
        return capillary_pressure_with_derivative().first;
    }
    // returns Pc and dPc/dln alpha_L
    [[nodiscard]] virtual std::pair<tmp<volScalarField>, tmp<volScalarField>>
    capillary_pressure_with_derivative() const = 0;

    inline static const word Pc_derivative_name = "dPc|dln(alpha_L)";

    // Declare runtime construction
    declareRunTimeSelectionTable(autoPtr, capillaryPressureModel, dictionary,
                                 (const dictionary &dict,
                                  const phaseInterface &interface),
                                 (dict, interface));
    // Selectors

    static autoPtr<capillaryPressureModel> New(const dictionary &dict,
                                               const phaseInterface &interface);
};

} // namespace Foam::capillaryPressureModels

#endif // MIKENO_multiphaseEuler_capillaryPressureModel_H
