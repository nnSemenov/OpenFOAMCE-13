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

Class
    Foam::heatTransferModels::Niegodajew

Description
    Niegodajew correlation for gas-liquid heat transfer in trickle bed. This
    correlation is specialized for tricked bed.

    Supports multiple partcicle phases, use sauter average.

    Reference:
    \verbatim
        P. Niegodajew, D. Asendrych / International Journal of Heat and Mass
        Transfer, 108(2017), 703-711.
    \endverbatim

Usage
    \table
        Property      | Description
        gas           | Gas phase name
        liquid        | Liquid phase name
        solids        | List of particle phases
        sigma         | Surface tension coefficient
        residualRe    | Lower bound for Re
    \endtable

    Example:
    \verbatim
    {
        type    Niegodajew;
        gas     air;
        liquid  water;
        solids  (solid);
        sigma   0.07;
    }
    \endverbatim

SourceFiles
    Gnielinski.C

\*---------------------------------------------------------------------------*/

#ifndef Niegodajew_H
#define Niegodajew_H

#include "heatTransferModel.H"
#include "averageDiameter.H"

namespace Foam::heatTransferModels
{

class Niegodajew : public heatTransferModel
{
  protected:
    //- Interface
    const dispersedPhaseInterface interface_;

    word gasName_;
    word liquidName_;
    List<word> solidNames_;

    scalar residualRe_;
    //- Surface tension
    dimensionedScalar sigma_;

  public:
    TypeName("Niegodajew");

    Niegodajew(const dictionary &dict, const phaseInterface &interface,
               const bool registerObject);

    [[nodiscard]] averageDiameter::sauterAverage average_diameter() const;

    [[nodiscard]] tmp<volScalarField> K(
        const scalar residualAlpha) const override;
};

} // namespace Foam::heatTransferModels

#endif // Niegodajew_H
