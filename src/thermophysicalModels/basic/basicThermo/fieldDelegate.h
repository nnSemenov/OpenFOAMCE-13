/*---------------------------------------------------------------------------*\
License
  This file is part of Mikeno. Mikeno (https://github.com/nnSemenov/Mikeno) is a
  fork of OpenFOAM (https://github.com/OpenFOAM/OpenFOAM-dev).

  Mikeno is free software: you can redistribute it and/or modify it  under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  Mikeno is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Mikeno.  If not, see <http://www.gnu.org/licenses/>.

Class
  Foam::fieldDelegate; Foam::volFieldDelegate

Description
  Field proxy classes to handle pressure offset

SourceFiles
  fieldDelegate.h

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_FIELDDELEGATE_H
#define MIKENO_FIELDDELEGATE_H

#include <cmath>
#include <vector>
#include <type_traits>

#include "volFields.H"
#include "physicalProperties.H"
#include "uniformGeometricFields.H"

namespace Foam
{

template <class F>
class fieldDelegate
{
  protected:
    const F &field;
    typename F::value_type offset;

  public:
    static_assert(std::is_same<typename F::value_type, scalar>::value,
                  "Offset must be scalar");

    fieldDelegate(const F &f, typename F::value_type of) : field(f), offset(of)
    {}

    auto operator[](label index) const { return this->field[index] + offset; }

    [[nodiscard]] label size() const { return this->field.size(); }
};

template <class vF>
class volFieldDelegate : public fieldDelegate<vF>
{
  public:
    using Patch = fieldDelegate<typename vF::Patch>;

  protected:
    std::vector<Patch> boundary;

  public:
    volFieldDelegate(const vF &f, typename vF::value_type of)
        : fieldDelegate<vF>(f, of)
    {
        const auto &boundaryField = this->field.boundaryField();
        forAll(boundaryField, patchIndex)
        {
            this->boundary.emplace_back(boundaryField[patchIndex],
                                        this->offset);
        }
    }

    const auto &boundaryField() const { return this->boundary; }
};

template <class vF>
inline auto absolutePressureVol(const vF &p, const dimensionedScalar &pOffset)
{
    return volFieldDelegate<vF>(p, pOffset.value());
}

template <class vF>
inline auto absolutePressureVol(const vF &p, const dimensionedScalar &pOffset,
                                const fvMesh &)
{
    return absolutePressureVol(p, pOffset);
}

inline auto absolutePressureVol(const uniformGeometricScalarField &p,
                                const dimensionedScalar &pOffset,
                                const fvMesh &mesh)
{
    //    using namespace Foam;
    const scalar value = p[0];
    if (not std::isfinite(p[0])) {
        return p;
    }
    return uniformGeometricScalarField(IOobject("pAbs", mesh.time().name(),
                                                mesh, IOobject::NO_READ,
                                                IOobject::NO_WRITE,
                                                false // don't register
                                                ),
                                       value + pOffset);
}

template <class pF>
inline auto absolutePressurePatch(const pF &p, const dimensionedScalar &pOffset)
{
    return fieldDelegate<pF>(p, pOffset.value());
}

template <class pF>
inline auto absolutePressurePatch(const pF &p, const dimensionedScalar &pOffset,
                                  const fvMesh &)
{
    return absolutePressurePatch(p, pOffset);
}

inline auto absolutePressurePatch(const UniformField<scalar> &p,
                                  const dimensionedScalar &pOffset,
                                  const fvMesh &)
{
    //    using namespace Foam;
    if (not std::isfinite(p[0])) {
        return p;
    }
    return UniformField<scalar>(p[0] + pOffset.value());
}

inline auto absolutePressurePatch(const UniformDimensionedField<scalar> &p,
                                  const dimensionedScalar &pOffset,
                                  const fvMesh &mesh)
{
    //    using namespace Foam;
    if (not std::isfinite(p[0])) {
        return p;
    }
    return UniformDimensionedField<scalar>(IOobject("pAbs", mesh.time().name(),
                                                    mesh, IOobject::NO_READ,
                                                    IOobject::NO_WRITE,
                                                    false // don't register
                                                    ),
                                           p[0] + pOffset);
}

} // namespace Foam

#endif // MIKENO_FIELDDELEGATE_H
