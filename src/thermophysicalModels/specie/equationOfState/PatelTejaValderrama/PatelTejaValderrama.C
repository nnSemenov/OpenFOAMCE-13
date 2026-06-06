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

#include "PatelTejaValderrama.H"
#include "IOstreams.H"
#include "cubicEOS.H"

// template <class Specie>
// Foam::PatelTejaValderrama<Specie>::PatelTejaValderrama(const word &name,
//                                                        const dictionary
//                                                        &dict)
//     : Specie{name, dict}, property_{dict}
//{
//     property_.requireRealGasEOS(name, true);
// }

template <class Specie>
void Foam::PatelTejaValderrama<Specie>::write(Ostream &os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    property_.write(dict);
    os << indent << dict.dictName() << dict;
}

template <class Specie>
Foam::Ostream &Foam::operator<<(Ostream &os,
                                const PatelTejaValderrama<Specie> &pg)
{
    pg.write(os);
    return os;
}

// template <class Specie>
// inline void Foam::PatelTejaValderrama<Specie>::operator*=(const scalar)
//{
//     noCoefficientMixing(PatelTejaValderrama);
// }
