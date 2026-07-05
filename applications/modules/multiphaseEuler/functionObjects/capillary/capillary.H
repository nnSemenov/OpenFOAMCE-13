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
    Foam::functionObjects::capillary

Description
    Writes capillary

Example
    \verbatim
    capillary{
        type capillary;
        libs ("libmultiphaseEulerFunctionObjects.so");

        writeCapillaryPressure  true;
        writePhasePressure      true;
        writeCapillaryForce     true;
    }
    \endvarbatim

SourceFiles
    capillary.C

\*---------------------------------------------------------------------------*/

#ifndef MIKENO_multiphaseEulerFunctionObjects_capillary_H
#define MIKENO_multiphaseEulerFunctionObjects_capillary_H

#include "fvMeshFunctionObject.H"
#include "phaseSystem.H"

#include "capillarySystem.H"

namespace Foam::functionObjects
{

class capillary : public fvMeshFunctionObject
{
  protected:
    bool write_capillary_pressure_;
    bool write_phase_pressure_;
    bool write_capillary_force_;

  public:
    TypeName("capillary");

    capillary(const word &name, const Time &runTime, const dictionary &dict);

    capillary(const capillary &) = delete;

    ~capillary() override = default;

    inline wordList fields() const override { return {}; }

    inline bool execute() override { return true; }

    bool write() override;

    bool read(const dictionary &dict) override;
};

} // namespace Foam::functionObjects

#endif // MIKENO_multiphaseEulerFunctionObjects_capillary_H
