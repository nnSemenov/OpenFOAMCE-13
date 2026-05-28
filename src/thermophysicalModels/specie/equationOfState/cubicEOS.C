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

#include "cubicEOS.H"
#include "dictionary.H"

#include <limits>

#include "../../../OpenFOAM/db/dictionary/dictionary.H"

using Foam::scalar;
using Foam::word;

Foam::realFluidProperty::realFluidProperty(const dictionary &dict)
{
    const dictionary &eos = dict.subDict("equationOfState");
    Tc_ = eos.lookup<scalar>("Tc");
    Vc_ = eos.lookup<scalar>("Vc");
    Pc_ = eos.lookup<scalar>("Pc");
    omega_ = eos.lookup<scalar>("omega");

    const word phase = eos.lookupOrDefault<word>("phase", "vapor");
    if (phase == "vapor") {
        is_gas = true;
    }
    else if (phase == "liquid") {
        is_gas = false;
    }
    else {
        FatalErrorInFunction << "Invalid value for phase: \"" << phase
                             << "\". Valid values: (vapor liquid)" << endl;
        ::Foam::abort(::Foam::FatalError);
    }
}

word Foam::realFluidProperty::checkForRealGasEOS(bool checkOmega) const noexcept
{

    if (this->Tc_ <= 0) {
        return "Invalid critical temperature: " + ::Foam::name(this->Tc_) +
            "[K]";
    }
    if (this->Pc_ <= 0) {
        return "Invalid critical pressure: " + ::Foam::name(this->Pc_) + "[Pa]";
    }
    if (this->Vc_ <= 0) {
        return "Invalid critical volume: " + ::Foam::name(this->Vc_) +
            "[m^3/kmol]";
    }
    if ((this->Zc() <= 0) or (this->Zc() > 1)) {
        return "Invalid critical compression factor: " +
            ::Foam::name(this->Pc_);
    }
    if (checkOmega and (this->omega_ <= -1)) {
        return "Invalid acentric factor: " + ::Foam::name(this->omega_);
    }
    return "";
}

void Foam::realFluidProperty::requireRealGasEOS(const word &specieName,
                                                bool require_omega) const
{
    const word error = checkForRealGasEOS(require_omega);
    if (error.empty()) {
        return;
    }

    FatalErrorInFunction << specieName << " is invalid specie: " << error
                         << endl;

    ::Foam::abort(::Foam::FatalError);
}

void Foam::realFluidProperty::write(dictionary &dict) const
{
    if (Tc_ > 0) {
        dict.add("Tc", Tc_);
    }
    if (Pc_ > 0) {
        dict.add("Pc", Pc_);
    }
    if (Vc_ > 0) {
        dict.add("Vc", Vc_);
    }
    if (omega_ > -1) {
        dict.add("omega", omega_);
    }
    const char *phase_val = this->is_gas ? "vapor" : "liquid";
    dict.add("phase", phase_val);
}

Foam::scalar Foam::solveCubicEquation(scalar a2, scalar a1, scalar a0,
                                      bool prefer_max, scalar lower_bound)
{
    const scalar Q = (3 * a1 - a2 * a2) / 9.0;
    const scalar Rl = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54.0;

    const scalar Q3 = Q * Q * Q;
    const scalar D = Q3 + Rl * Rl;

    //  scalar root = -1;

    if (D <= 0) {
        // Three roots are real, but may contains negative
        const scalar th = ::acos(Rl / sqrt(-Q3));
        const scalar qm = 2 * sqrt(-Q);
        const scalar r1 = qm * cos(th / 3.0) - a2 / 3.0;
        const scalar r2 =
            qm * cos((th + 2 * constant::mathematical::pi) / 3.0) - a2 / 3.0;
        const scalar r3 =
            qm * cos((th + 4 * constant::mathematical::pi) / 3.0) - a2 / 3.0;

        if (prefer_max) {
            return max(r1, max(r2, r3));
        }
        return take_min_valid(r1, r2, r3, lower_bound);
    }
    else {
        // One root is real
        const scalar D05 = sqrt(D);
        const scalar S = cbrt(Rl + D05);
        scalar Tl = 0;
        if (D05 > Rl) {
            Tl = -cbrt(mag(Rl - D05));
        }
        else {
            Tl = cbrt(Rl - D05);
        }

        const scalar root = S + Tl - a2 / 3.0;

        return root;
    }
}

word Foam::parseBinarySpeciePair(const word &str, word &sp0, word &sp1)
{
    sp0 = "";
    sp1 = "";
    const word explain =
        "it is not a valid binary specie pair. Example: \"O2:CO2\"";
    //    const std::string_view str{str__};
    const size_t comma_idx = str.find_first_of(':');
    if (comma_idx == str.npos) {
        return str + " doesn't contain colon, " + explain;
    }
    if (comma_idx not_eq str.find_last_of(':')) {
        return str + " has multiple colon, " + explain;
    }
    sp0 = std::string{str.data(), str.data() + comma_idx};
    sp1 = std::string{str.begin() + comma_idx + 1, str.end()};
    return "";
}

void Foam::parseBinaryInteractionMatrix(
    const dictionary &kDict, const std::function<bool(const word &)> &skipKey,
    scalarSquareMatrix &mat, const speciesTable &spTable, bool symmetry)
{
    for (label r = 0; r < mat.m(); r++) {
        for (label c = 0; c < mat.n(); c++) {
            mat(r, c) = std::numeric_limits<scalar>::signaling_NaN();
        }
    }

    auto set_k_ij = [&mat](label i, label j, scalar k) noexcept {
        const bool alreadySet = std::isfinite(mat(i, j));
        if (alreadySet) {
            WarningInFunction << "Duplicated declaration of k(" << i << ',' << j
                              << "), previous value will be covered." << endl;
        }
        mat(i, j) = k;
    };
    const auto keys = kDict.keys(true);
    forAll(keys, keyIdx)
    {
        const word &speciePair = keys[keyIdx];
        if (skipKey(speciePair)) {
            continue;
        }
        word sp0, sp1;
        const word err = parseBinarySpeciePair(speciePair, sp0, sp1);
        if (not err.empty()) {
            FatalErrorInFunction << "Invalid binary interaction table: " << err
                                 << endl;
            ::Foam::abort(::Foam::FatalError);
        }
        const label idxSp0 = spTable[sp0];
        const label idxSp1 = spTable[sp1];

        const scalar k = kDict.template lookup<scalar>(speciePair);
        set_k_ij(idxSp0, idxSp1, k);
        if (symmetry) {
            set_k_ij(idxSp1, idxSp0, k);
        }
    }

    for (label r = 0; r < mat.m(); r++) {
        for (label c = 0; c < mat.n(); c++) {
            scalar &k = mat(r, c);
            if (not std::isfinite(k)) {
                k = 0;
            }
        }
    }
}