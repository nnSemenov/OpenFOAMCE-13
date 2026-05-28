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

#ifndef effectiveHeatConductivity_H
#define effectiveHeatConductivity_H

#include <span>
#include "volFieldsFwd.H"
#include "typeInfo.H"

namespace Foam
{

class effectiveHeatConductivity
{
  protected:
    scalar alpha_small_{1e-6};

  public:
    // TypeName("effectiveHeatConducitvity");

    explicit effectiveHeatConductivity(const dictionary &dict) {}

    virtual ~effectiveHeatConductivity() = default;

    struct kappaInfo
    {
        const volScalarField &alpha;
        tmp<volScalarField> kappa;
    };

    virtual tmp<volScalarField> kappaEff(std::span<kappaInfo>) const = 0;

    scalar alpha_small() const { return this->alpha_small_; }

    scalar &alpha_small() { return this->alpha_small_; }

    static autoPtr<effectiveHeatConductivity> New(const dictionary &dict);
};

} // namespace Foam

#endif // effectiveHeatConductivity_H
