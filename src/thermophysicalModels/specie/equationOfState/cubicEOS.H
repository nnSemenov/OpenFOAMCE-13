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
    Foam::realFluidProperty

Description
    Utilities for real fluid EOS.

SourceFiles
    cubicEOS.C
\*---------------------------------------------------------------------------*/

#ifndef CUBICEOS_H
#define CUBICEOS_H

#include <cassert>
#include <type_traits>
#include <functional>

#include "autoPtr.H"
#include "mathematicalConstants.H"
#include "scalarMatrices.H"
#include "speciesTable.H"
#include "thermodynamicConstants.H"

namespace Foam
{
using constant::thermodynamic::RR;

struct realFluidProperty
{

    //- Critical Temperature [K]
    scalar Tc_{-1};

    //- Critical volume [m^3/kmol]
    scalar Vc_{-1};

    //- Critical Pressure [Pa]
    scalar Pc_{-1};

    //- Acentric factor [-]
    scalar omega_{-2};

    bool is_gas{true};

    realFluidProperty() = default;

    explicit realFluidProperty(const dictionary &dict);

    [[nodiscard]] inline scalar Zc() const noexcept
    {
        return Pc_ * Vc_ / (RR * Tc_);
    }

    [[nodiscard]] word checkForRealGasEOS(bool require_omega) const noexcept;

    void requireRealGasEOS(const word &specieName, bool require_omega) const;

    void write(dictionary &eos_dict) const;
};

struct cubicEOSDimlessCoefficient
{
    scalar A;
    scalar B;
};

// Solve z^3 + a2 *z^2 + a1 * z + a0 = 0
scalar solveCubicEquation(scalar a2, scalar a1, scalar a0, bool prefer_max,
                          scalar lower_bound);

constexpr inline scalar take_min_valid(scalar x0, scalar x1, scalar x2,
                                       scalar lower_bound) noexcept
{
    constexpr scalar gr = std::numeric_limits<scalar>::max();

    const scalar r0 = (x0 > lower_bound ? x0 : gr);
    const scalar r1 = (x1 > lower_bound ? x1 : gr);
    const scalar r2 = (x2 > lower_bound ? x2 : gr);
    const scalar r = std::min(r0, std::min(r1, r2));
    return r;
}

struct residualEnergy
{
    scalar h_R; //[J/kg]
    scalar e_R; //[J/kg]
    scalar Z;
};

struct residualSpecificHeat
{
    scalar Cp_R;            //[J/kg/K]
    scalar Cv_R;            //[J/kg/K]
    scalar Cp_R_minus_Cv_R; //[J/kg/K]
    scalar Z;
};

struct residualSpecificHeatDerivative
{
    scalar dCpdT_R;           //[J/kg/K^2]
    scalar dCvdT_R;           //[J/kg/K^2]
    scalar dCp_minus_Cv_dT_R; //[J/kg/K^2]
    scalar Z;
};

namespace internal
{
// Some SFINAE tricks
template <class EOS>
inline int requiresCubicEOS(const EOS *eos)
{
    constexpr scalar T = 273;
    constexpr scalar p = 1e5;
    constexpr scalar W = 28.0;
    // EOS must have these interfaces
    [[maybe_unused]] typename EOS::EOSCore core = eos->core(p, T);
    [[maybe_unused]] const scalar u = core.Z(p, T) + core.rho(p, T, W) +
        core.h(p, T, W) + core.Cp(p, T, W) + core.e(p, T, W) +
        core.Cv(p, T, W) + core.sp(p, T, W) + core.sv(p, T, W) +
        core.psi(p, T, W) + core.CpCv(p, T, W) + core.alphav(p, T, W);
    return 0;
}

inline void requiresCubicEOS(...) {}

} // namespace internal
template <class EOS>
struct is_cubic_EOS
{
    static constexpr bool value =
        std::is_same_v<decltype(internal::requiresCubicEOS(
                           static_cast<const EOS *>(nullptr))),
                       int> and
        (not EOS::incompressible);
};

word parseBinarySpeciePair(const word &str, word &sp0, word &sp1);

void parseBinaryInteractionMatrix(
    const dictionary &coeffDict,
    const std::function<bool(const word &)> &skipKey, scalarSquareMatrix &mat,
    const speciesTable &spTable, bool symmetry);
} // namespace Foam

#endif // CUBICEOS_H
