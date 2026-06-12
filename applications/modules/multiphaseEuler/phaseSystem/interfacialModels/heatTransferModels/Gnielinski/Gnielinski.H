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
    Foam::heatTransferModels::Gnielinski

Description
    Gnielinski (1975~1982) correlation for packed bed heat transfer. Suitable
    for both laminar and turbulence. If spherical, form factor is related to
    porosity (this correlation has better handling for identical spheres);
    otherwise form factor is input constant.

    Factors for non-spherical particles: 1.6 for cylindrical particles with
    0.24<l/d<1.2; 1.6 for cubes; 2.1 for Raschig rings; 2.3 for Berl saddles.

    factor = 1+1.5*(1-porosity) for spherical particles.

    Reference:
    \verbatim

        Gnielinski, V. (1975). Berechnung mittlerer Wärme- und
        Stoffübergangskoeffizienten an laminar und turbulent überströmten
        Einzelkörpern mit Hilfe einer einheitlichen Gleichung. Forsch Ing-Wes,
        41, 145–153.

        Gnielinski, V. (1981). Equations for the calculation of heat and mass
        transfer during flow through stationary spherical packings at moderate
        and high Peclet numbers. Int Chem Eng, 21(3), 378–383.

        Gnielinski, V. (1982). Berechnung des Wärme- und Stoffaustauschs in
        durchströmten ruhenden Schüttungen. Verfahrenstechnik, 16(1), 36–39.

        VDI e. V. (2010). VDI Heat Atlas. Springer Berlin, Heidelberg.
        https://doi.org/10.1007/978-3-540-77877-6
    \endverbatim

Usage
    \table
        Property      | Description
        sphere        | Whether particles are sphere.
        factor        | Form factor for non-spherical particles.
        residualRe    | Lower bound for Re. Defaults to 1e-6.
    \endtable

    Example:
    \verbatim
    {
        type    Gnielinski;
        sphere  yes;
    }
    \endverbatim

SourceFiles
    Gnielinski.C

\*---------------------------------------------------------------------------*/

#ifndef Gnielinski_H
#define Gnielinski_H

#include "heatTransferModel.H"

namespace Foam::heatTransferModels
{

class Gnielinski : public heatTransferModel
{
  public:

  protected:
    //- Interface
    const dispersedPhaseInterface interface_;
    // True if is sphere
    bool sphere_;
    // Shape factor for non-spherical particles. Not used when sphere
    scalar shapeFactor_;
    // Lower bound of Re.
    scalar residualRe_{1e-6};

  public:
    TypeName("Gnielinski");

    Gnielinski(const dictionary &dict, const phaseInterface &interface,
               const bool registerObject);

    virtual ~Gnielinski() override = default;

    [[nodiscard]] inline scalar residualRe() const noexcept
    {
        return residualRe_;
    }
    //- The heat transfer function K used in the enthalpy equation
    tmp<volScalarField> K(const scalar residualAlpha) const override;
};

} // namespace Foam::heatTransferModels

#endif // Gnielinski_H
