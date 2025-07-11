/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::hexRef8

Description
    Refinement of (split) hexes using polyTopoChange.

SourceFiles
    hexRef8.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef8_H
#define hexRef8_H

#include "labelIOList.H"
#include "face.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "primitivePatch.H"
#include "removeFaces.H"
#include "refinementHistory.H"
#include "PackedBoolList.H"
#include "uniformDimensionedFields.H"
#include "cellShapeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;
class polyTopoChange;
class polyTopoChangeMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                           Class hexRef8 Declaration
\*---------------------------------------------------------------------------*/

class hexRef8
{
    // Private Data

        //- Reference to underlying mesh.
        const polyMesh& mesh_;

        //- Per cell the refinement level
        labelIOList cellLevel_;

        //- Per point the refinement level
        labelIOList pointLevel_;

        //- Typical edge length between unrefined points
        uniformDimensionedScalarField level0Edge_;

        //- Refinement history
        refinementHistory history_;

        //- Face remover engine
        removeFaces faceRemover_;

        //- Level of saved points
        Map<label> savedPointLevel_;

        //- Level of saved cells
        Map<label> savedCellLevel_;

        //- Cell shapes when seen as split hexes
        mutable autoPtr<cellShapeList> cellShapesPtr_;


    // Private Member Functions

        //- Reorder according to map.
        static void reorder
        (
            const labelList& map,
            const label len,
            const label null,
            labelList& elems
        );

        //- Get patch index
        label getPatchIndex(const label facei) const;

        //- Adds a face on top of existing facei. Reverses if necessary.
        label addFace
        (
            polyTopoChange& meshMod,
            const label facei,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Adds internal face from point. No checks on reversal.
        //  Adds an internal face from an edge. Assumes orientation correct.
        //  Problem is that the face is between four new vertices. So what do
        //  we provide as master? The only existing mesh item we have is the
        //  edge we have split.  Have to be careful in only using it if it has
        //  internal faces since otherwise polyMeshMorph will complain (because
        //  it cannot generate a sensible mapping for the face)
        label addInternalFace
        (
            polyTopoChange& meshMod,
            const label meshFacei,
            const label meshPointi,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Modifies existing facei for either new owner/neighbour or new face
        //  points. Reverses if necessary.
        void modifyFace
        (
            polyTopoChange& meshMod,
            const label facei,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Bit complex way to determine the unrefined edge length.
        scalar getLevel0EdgeLength() const;

        //- Get cell added to point of celli (if any)
        //  Check whether pointi is an anchor on celli.  If it is not check
        //  whether any other point on the face is an anchor cell.
        label getAnchorCell
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label celli,
            const label facei,
            const label pointi
        ) const;

        //- Get new owner and neighbour (in unspecified order) of pointi
        //  on facei.
        void getFaceNeighbours
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label facei,
            const label pointi,

            label& own,
            label& nei
        ) const;


        //- Get index of minimum pointlevel.
        label findMinLevel(const labelList& f) const;

        //- Get maximum pointlevel.
        label findMaxLevel(const labelList& f) const;

        //- Count number of vertices <= anchorLevel
        label countAnchors(const labelList&, const label) const;

        //- Debugging: dump cell as .obj file
        void dumpCell(const label celli) const;

        //- Find index of point with wantedLevel, starting from fp.
        label findLevel
        (
            const label facei,
            const face& f,
            const label startFp,
            const bool searchForward,
            const label wantedLevel
        ) const;

        //- debug:check orientation of added internal face
        static void checkInternalOrientation
        (
            polyTopoChange& meshMod,
            const label celli,
            const label facei,
            const point& ownPt,
            const point& neiPt,
            const face& newFace
        );

        //- debug:check orientation of new boundary face
        static void checkBoundaryOrientation
        (
            polyTopoChange& meshMod,
            const label celli,
            const label facei,
            const point& ownPt,
            const point& boundaryPt,
            const face& newFace
        );

        //- If p0 and p1 are existing vertices check if edge is split and insert
        // splitPoint.
        void insertEdgeSplit
        (
            const labelList& edgeMidPoint,
            const label p0,
            const label p1,
            DynamicList<label>& verts
        ) const;

        //- Store in maps correspondence from midpoint to anchors and faces.
        //  Internal faces are one per edge between anchor points. So one per
        //  midPoint between the anchor points. Here we store the information on
        //  the midPoint and if we have enough information:
        //  - two anchors
        //  - two face mid points
        //  we add the face. Note that this routine can get called anywhere from
        //  two times (two unrefined faces) to four times (two refined faces) so
        //  the first call that adds the information creates the face.
        label storeMidPointInfo
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& edgeMidPoint,
            const label celli,
            const label facei,
            const bool faceOrder,
            const label midPointi,
            const label anchorPointi,
            const label faceMidPointi,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod
        ) const;

        //- Create all internal faces from an unsplit face.
        void createInternalFromSplitFace
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& edgeMidPoint,
            const label celli,
            const label facei,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod,
            label& nFacesAdded
        ) const;

        //- Create all internal faces to split celli into 8.
        void createInternalFaces
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& faceAnchorLevel,
            const labelList& edgeMidPoint,
            const label celli,
            polyTopoChange& meshMod
        ) const;

        //- Store vertices from startFp up to face split point.
        //  Used when splitting face into 4.
        void walkFaceToMid
        (
            const labelList& edgeMidPoint,
            const label cLevel,
            const label facei,
            const label startFp,
            DynamicList<label>& faceVerts
        ) const;

        //- Same as walkFaceToMid but now walk back.
        void walkFaceFromMid
        (
            const labelList& edgeMidPoint,
            const label cLevel,
            const label facei,
            const label startFp,
            DynamicList<label>& faceVerts
        ) const;

        //- Updates refineCell so consistent 2:1 refinement. Returns local
        //  number of cells changed.
        label faceConsistentRefinement
        (
            const bool maxSet,
            PackedBoolList& refineCell
        ) const;

        //- Check wanted refinement for 2:1 consistency
        void checkWantedRefinementLevels(const labelList&) const;


        // Cellshape recognition

            //- Collect all points on face of certain level
            void collectLevelPoints
            (
                const labelList& f,
                const label level,
                DynamicList<label>& points
            ) const;

            //- Collect all points on face (in local numbering) of certain level
            void collectLevelPoints
            (
                const labelList& meshPoints,
                const labelList& f,
                const label level,
                DynamicList<label>& points
            ) const;

            //- Collect all faces with four corner points and return true if
            //  hex was matched (6 faces of each four corner points)
            bool matchHexShape
            (
                const label celli,
                const label cellLevel,
                DynamicList<face>& quads
            ) const;


public:

    //- Runtime type information
    ClassName("hexRef8");


    // Constructors

        //- Construct from mesh, read_if_present refinement data
        //  (from write below). If readHistory is true does read_if_present
        //  of refinement history. If false clears all history
        hexRef8(const polyMesh& mesh, const bool readHistory = true);

        //- Construct from mesh and un/refinement data and optional size of
        //  starting cells
        hexRef8
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const refinementHistory& history,
            const scalar level0Edge = -great
        );

        //- Construct from mesh and refinement data and optional size of
        //  starting cells
        hexRef8
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar level0Edge = -great
        );

        //- Disallow default bitwise copy construction
        hexRef8(const hexRef8&) = delete;


    // Member Functions

        // Access

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            const labelIOList& cellLevel() const
            {
                return cellLevel_;
            }

            const labelIOList& pointLevel() const
            {
                return pointLevel_;
            }

            const refinementHistory& history() const
            {
                return history_;
            }

            //- Typical edge length between unrefined points
            scalar level0EdgeLength() const
            {
                return level0Edge_.value();
            }

        // Refinement

            //- Gets level such that the face has four points <= level.
            label faceLevel(const label facei) const;

            //- Given valid mesh and current cell level and proposed
            //  cells to refine calculate any clashes (due to 2:1) and return
            //  ok list of cells to refine.
            //  Either adds cells to refine to set (maxSet = true) or
            //  removes cells to refine (maxSet = false)
            labelList consistentRefinement
            (
                const labelList& cellsToRefine,
                const bool maxSet
            ) const;

            //- Like consistentRefinement but slower:
            //
            //  - specify number of cells between consecutive refinement levels
            //    (consistentRefinement equivalent to 1)
            //  - specify max level difference between point-connected cells.
            //    (-1 to disable) Note that with normal 2:1 limitation
            //    (maxFaceDiff=1) there can be 8:1 size difference across point
            //    connected cells so maxPointDiff allows you to make that less.
            //  cellsToRefine : cells we're thinking about refining. It will
            //                  extend this set. All refinement levels will be
            //                  at least maxFaceDiff layers thick.
            //  facesToCheck  : additional faces where to implement the
            //                  maxFaceDiff thickness (usually only boundary
            //                  faces)
            labelList consistentSlowRefinement
            (
                const label maxFaceDiff,
                const labelList& cellsToRefine,
                const labelList& facesToCheck,
                const label maxPointDiff,
                const labelList& pointsToCheck
            ) const;

            //- Like consistentSlowRefinement but uses different meshWave
            //  (proper distance instead of topological count). No point checks
            //  yet.
            labelList consistentSlowRefinement2
            (
                const label maxFaceDiff,
                const labelList& cellsToRefine,
                const labelList& facesToCheck
            ) const;

            //- Insert refinement. All selected cells will be split into 8.
            //  Returns per element in cells the 8 cells they were split into.
            //  Guarantees that the 0th element is the original cell label.
            //  Mapping:
            //  -split cells: 7 new ones get added from original
            //  -split faces: original gets modified; 3 new ones get added
            //               from original
            //  -added internal faces: added from original cell face(if
            //   that was internal) or created out-of-nothing (so will not
            //   get mapped!).
            //  -points added to split edge: added from edge start()
            //  -midpoints added: added from cellPoints[0].
            labelListList setRefinement
            (
                const labelList& cells,
                polyTopoChange&
            );

            //- Update local numbering for changed mesh.
            //  Called after the mesh change. setRefinement will already have
            //  made sure the pointLevel_ and cellLevel_ are the size of the new
            //  mesh so we only need to account for reordering.
            void topoChange(const polyTopoChangeMap&);


            // Restoring : is where other processes delete and reinsert data.
            //             These callbacks allow this to restore the cellLevel
            //             and pointLevel for reintroduced points.
            //             Is not related to undoing my refinement

                //- Signal points/face/cells for which to store data
                void storeData
                (
                    const labelList& pointsToStore,
                    const labelList& facesToStore,
                    const labelList& cellsToStore
                );

                //- Update local numbering + undo
                //  Data to restore given as new pointlabel + stored pointlabel
                //  (i.e. what was in pointsToStore)
                void topoChange
                (
                    const polyTopoChangeMap&,
                    const Map<label>& pointsToRestore,
                    const Map<label>& facesToRestore,
                    const Map<label>& cellsToRestore
                );


            //- Update local numbering for subsetted mesh.
            //  Gets new-to-old maps. Not compatible with unrefinement.
            void subset
            (
                const labelList& pointMap,
                const labelList& faceMap,
                const labelList& cellMap
            );

            //- Update local numbering for mesh redistribution
            void distribute(const polyDistributionMap&);

            //- Debug: Check coupled mesh for correctness
            void checkMesh() const;

            //- Debug: Check 2:1 consistency across faces.
            //  maxPointDiff==-1 : only check 2:1 across faces
            //  maxPointDiff!=-1 : check point-connected cells.
            void checkRefinementLevels
            (
                const label maxPointDiff,
                const labelList& pointsToCheck
            ) const;

            //- Utility: get hexes as cell shapes
            const cellShapeList& cellShapes() const;


        // Unrefinement (undoing refinement, not arbitrary coarsening)

            //- Return the points at the centre of top-level split cells
            //  that can be unsplit.
            labelList getSplitPoints() const;

            //- Given proposed
            //  splitPoints to unrefine according to calculate any clashes
            //  (due to 2:1) and return ok list of points to unrefine.
            //  Either adds points to refine to set (maxSet = true) or
            //  removes points to refine (maxSet = false)
            labelList consistentUnrefinement
            (
                const labelList& pointsToUnrefine,
                const bool maxSet
            ) const;

            //- Remove some refinement. Needs to be supplied output of
            //  consistentUnrefinement. Only call if undoable set.
            //  All 8 pointCells of a split point will be combined into
            //  the lowest numbered cell of those 8.
            void setUnrefinement
            (
                const labelList& splitPointLabels,
                polyTopoChange&
            );

        // Write

            // Set instance for mesh files
            void setInstance(const fileName& inst);

            //- Force writing refinement+history to polyMesh directory.
            bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hexRef8&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
