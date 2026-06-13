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

Namespace
    Foam

Description
    Average diameter of multiple particle phases.
\*---------------------------------------------------------------------------*/

#ifndef averageDiameter_H
#define averageDiameter_H

#include <cassert>
#include "phaseModel.H"
#include "volFields.H"

namespace Foam::averageDiameter
{
struct sauterAverage
{
    //- Sum of particle volume fraction
    volScalarField alpha;
    //- Average diameter
    volScalarField d;
};

template <class getPhaseModelFun>
[[nodiscard]] inline sauterAverage sauter(label nSolidPhases,
                                          getPhaseModelFun getPhaseModel)
{
    assert(nSolidPhases > 0);
    const phaseModel &phase0 = getPhaseModel(0);
    volScalarField alpha_sum(max(static_cast<const volScalarField &>(phase0),
                                 phase0.residualAlpha()));
    volScalarField inv_diameter_sum(alpha_sum / phase0.d());

    for (label i = 1; i < nSolidPhases; i++) {
        const phaseModel &phase = getPhaseModel(i);
        volScalarField alpha_i(max(static_cast<const volScalarField &>(phase),
                                   phase.residualAlpha()));
        alpha_sum += alpha_i;
        inv_diameter_sum += alpha_i / phase.d();
    }
    inv_diameter_sum = alpha_sum / inv_diameter_sum;
    return {
        .alpha = std::move(alpha_sum),
        .d = std::move(inv_diameter_sum),
    };
}

} // namespace Foam::averageDiameter

#endif // averageDiameter_H
