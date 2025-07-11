/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::heatTransfer

Description
    Model for heat exchange. Requires specification of an ambient temperature
    with which to exchange heat, and a model for the heat transfer coefficient
    (htc) and the area per unit volume (Av). These are then used to apply the
    following source to the energy equation:

    \f[
        -htc*Av*(T_a - T)
    \f]

    If the semiImplicit option is set, then this becomes:

    \f[
        -htc*Av*(T_a - T) + htc*Av/Cp*h - Sp(htc*Av/Cp, h);
    \f]

Usage
    Example usage:
    \verbatim
    heatTransfer
    {
        type            heatTransfer;

        cellZone        c0;

        semiImplicit    no;

        Ta              300;

        Av              200;

        heatTransferCoefficientModel constant;

        htc             10;
    }
    \endverbatim

See also
    Foam::fv::heatTransferCoefficientModel

SourceFiles
    heatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransfer_H
#define heatTransfer_H

#include "fvModel.H"
#include "fvCellZone.H"
#include "heatTransferCoefficientModel.H"
#include "heatTransferAv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class heatTransfer Declaration
\*---------------------------------------------------------------------------*/

class heatTransfer
:
    public fvModel
{
    // Private data

        //- The cellzone the model applies to
        fvCellZone zone_;

        //- Name of the phase
        word phaseName_;

        //- Flag to activate semi-implicit coupling
        bool semiImplicit_;

        //- Name of temperature field; default = "T"
        word TName_;

        //- Ambient temperature
        dimensionedScalar Ta_;

        //- The heat transfer area per unit volume
        autoPtr<heatTransferAv> heatTransferAv_;

        //- The heat transfer model
        autoPtr<heatTransferCoefficientModel> heatTransferCoefficientModel_;


    // Private member functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Source term to energy equation
        template<class AlphaFieldType>
        inline void add
        (
            const AlphaFieldType& alpha,
            fvMatrix<scalar>& eqn
        ) const;


public:

    //- Runtime type information
    TypeName("heatTransfer");


    // Constructors

        //- Construct from dictionary
        heatTransfer
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~heatTransfer();


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Source term to energy equation
            virtual void addSup
            (
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Source term to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;


        // Correction

            //- Correct the model
            virtual void correct();


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
