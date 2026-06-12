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
    Foam::heatTransferModels::PowerLawNu

Description
    Traditional power-law heat transfer correlation. Nu = a * Re^b * Pr^c.

Usage
    \table
        Property      | Description
        a             | Coefficient to Nu [-]
        b             | Exponent to Re [-]
        c             | Exponent to Pr [-]
        superficialRe | Use superficial velocity for Re. Defaults to false.
        residualRe    | Lower bound of Re for stabilization. Defaults to 1e-6;
    \endtable

    Example:
    \verbatim
    {
        type PowerLawNu;
        a 0.023;
        b 0.8;
        c 0.4;
        residualRe 1e-6;
    }
    \endverbatim

SourceFiles
    PowerLawNu.C

\*---------------------------------------------------------------------------*/

#ifndef PowerLawNu_H
#define PowerLawNu_H

#include "heatTransferModel.H"

namespace Foam::heatTransferModels
{

class PowerLawNu : public heatTransferModel
{

  protected:
    //- Interface
    const dispersedPhaseInterface interface_;

    dimensionedScalar coeff_a_;
    dimensionedScalar coeff_b_;
    dimensionedScalar coeff_c_;
    bool superficial_Re_{false};
    scalar residual_Re_{1e-6};

  public:
    TypeName("PowerLawNu");

    PowerLawNu(const dictionary &dict, const phaseInterface &interface,
               const bool registerObject);

    [[nodiscard]] inline auto &coeff_a() const noexcept { return coeff_a_; }
    [[nodiscard]] inline auto &coeff_b() const noexcept { return coeff_b_; }
    [[nodiscard]] inline auto &coeff_c() const noexcept { return coeff_c_; }
    [[nodiscard]] inline bool superficialRe() const noexcept
    {
        return superficial_Re_;
    }

    [[nodiscard]] inline scalar residualRe() const noexcept
    {
        return residual_Re_;
    }

    virtual ~PowerLawNu() override = default;

    //- The heat transfer function K used in the enthalpy equation
    tmp<volScalarField> K(const scalar residualAlpha) const override;
};
} // namespace Foam::heatTransferModels

#endif