lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(em edgeMesh)
add_library(meshTools SHARED
    ${em}/edgeMesh.C
    ${em}/edgeMeshIO.C
    ${em}/edgeMeshNew.C
)

set(edgeMeshFormats ${em}/edgeMeshFormats)
target_sources(meshTools PRIVATE
    ${edgeMeshFormats}/edgeMeshFormatsCore.C

    ${edgeMeshFormats}/edgeMesh/edgeMeshFormat.C
    ${edgeMeshFormats}/edgeMesh/edgeMeshFormatRunTime.C

    ${edgeMeshFormats}/extendedFeatureEdgeMesh/extendedFeatureEdgeMeshFormat.C
    ${edgeMeshFormats}/extendedFeatureEdgeMesh/extendedFeatureEdgeMeshFormatRunTime.C

    ${edgeMeshFormats}/nas/NASedgeFormat.C
    ${edgeMeshFormats}/nas/NASedgeFormatRunTime.C

    ${edgeMeshFormats}/obj/OBJedgeFormat.C
    ${edgeMeshFormats}/obj/OBJedgeFormatRunTime.C

    ${edgeMeshFormats}/starcd/STARCDedgeFormat.C
    ${edgeMeshFormats}/starcd/STARCDedgeFormatRunTime.C

    ${edgeMeshFormats}/vtk/VTKedgeFormat.C
    ${edgeMeshFormats}/vtk/VTKedgeFormatRunTime.C

    ${em}/featureEdgeMesh/featureEdgeMesh.C
)

set(eem ${em}/extendedEdgeMesh)
target_sources(meshTools PRIVATE

    ${eem}/extendedEdgeMesh.C
    ${eem}/extendedEdgeMeshNew.C

    ${eem}/extendedEdgeMeshFormats/extendedEdgeMeshFormat/extendedEdgeMeshFormat.C
    ${eem}/extendedEdgeMeshFormats/extendedEdgeMeshFormat/extendedEdgeMeshFormatRunTime.C

)

set(efm ${eem}/extendedFeatureEdgeMesh)
target_sources(meshTools PRIVATE
    ${efm}/extendedFeatureEdgeMesh.C
)

target_sources(meshTools PRIVATE
    cellClassification/cellClassification.C
    cellClassification/cellInfo.C

    patchDist/patchDistWave/patchDistWave.C

    cellFeatures/cellFeatures.C
)

set(csys coordinateSystems)
target_sources(meshTools PRIVATE
    ${csys}/coordinateSystem.C
    ${csys}/coordinateSystemNew.C
    ${csys}/coordinateSystems.C
    ${csys}/cylindricalCS.C
    ${csys}/cartesianCS.C
    ${csys}/coordinateRotation/axesRotation.C
    ${csys}/coordinateRotation/coordinateRotation.C
    ${csys}/coordinateRotation/coordinateRotationNew.C
    ${csys}/coordinateRotation/EulerCoordinateRotation.C
    ${csys}/coordinateRotation/STARCDCoordinateRotation.C
    ${csys}/coordinateRotation/cylindrical.C
)

target_sources(meshTools PRIVATE
    edgeFaceCirculator/edgeFaceCirculator.C

    polyMeshZipUpCells/polyMeshZipUpCells.C

    meshSearch/meshSearchBoundBox.C
    meshSearch/meshSearch.C
    meshSearch/meshBoundarySearch.C

    meshTools/meshTools.C

    algorithms/FaceCellWave/FaceCellWaveName.C
    algorithms/PointEdgeWave/PointEdgeWaveName.C
    algorithms/PatchEdgeFaceWave/PatchEdgeFaceWaveName.C
    algorithms/PatchEdgeFaceWave/patchPatchDist.C

    regionSplit/regionSplit.C
    regionSplit/localPointRegion.C

    indexedOctree/treeDataEdge.C
    indexedOctree/treeDataFace.C
    indexedOctree/treeDataPoint.C
    indexedOctree/treeDataPrimitivePatchName.C
    indexedOctree/treeDataTriSurface.C

    searchableSurfaces/searchableSurface/searchableSurface.C
    searchableSurfaces/searchableSurfaceList/searchableSurfaceList.C

    searchableSurfaces/box/box_searchableSurface.C
    searchableSurfaces/cylinder/cylinder_searchableSurface.C
    searchableSurfaces/disk/disk_searchableSurface.C
    searchableSurfaces/plane/plane_searchableSurface.C
    searchableSurfaces/plate/plate_searchableSurface.C
    searchableSurfaces/sphere/sphere_searchableSurface.C
    searchableSurfaces/withGaps/withGaps_searchableSurface.C
    searchableSurfaces/triSurface/triSurface_searchableSurface.C
    searchableSurfaces/closedTriSurface/closedTriSurface.C
    searchableSurfaces/extrudedCircle/extrudedCircle_searchableSurface.C
    searchableSurfaces/collection/collection_searchableSurface.C

    searchableSurfaces/searchableSurfacesQueries/searchableSurfacesQueries.C
    searchableSurfaces/triSurface/extractCloseness.C

    topoSets/cellSet.C
    topoSets/topoSet.C
    topoSets/faceSet.C
    topoSets/pointSet.C

    momentOfInertia/momentOfInertia.C

    triSurface/orientedSurface/orientedSurface.C
    triSurface/surfaceLocation/surfaceLocation.C
)

set(booleanOps triSurface/booleanOps)
set(surfaceIntersection ${booleanOps}/surfaceIntersection)
target_sources(meshTools PRIVATE
    ${surfaceIntersection}/surfaceIntersection.C
    ${surfaceIntersection}/surfaceIntersectionFuncs.C
    ${surfaceIntersection}/edgeIntersections.C
)

set(booleanSurface  ${booleanOps}/booleanSurface)
target_sources(meshTools PRIVATE
    ${booleanSurface}/booleanSurface.C
)

set(intersectedSurface  ${booleanOps}/intersectedSurface)
target_sources(meshTools PRIVATE

    ${intersectedSurface}/intersectedSurface.C
    ${intersectedSurface}/edgeSurface.C
)

target_sources(meshTools PRIVATE
    triSurface/triSurfaceSearch/triSurfaceSearch.C
    triSurface/triSurfaceSearch/triSurfaceRegionSearch.C
    triSurface/triangleFuncs/triangleFuncs.C
    triSurface/surfaceFeatures/surfaceFeatures.C
    triSurface/triSurfaceTools/triSurfaceTools.C
    triSurface/triSurfaceTools/geompack/geompack.C
    triSurface/triSurfaceTools/pointToPointPlanarInterpolation.C

    twoDPointCorrector/twoDPointCorrector.C

    patchToPatchTools/patchToPatchTools.C

    patchToPatch/patchToPatch/patchToPatch.C
    patchToPatch/patchToPatch/patchToPatchParallelOps.C
    patchToPatch/matching/matchingPatchToPatch.C
    patchToPatch/nearby/nearbyPatchToPatch.C
    patchToPatch/nearest/nearestPatchToPatch.C
    patchToPatch/inverseDistance/inverseDistancePatchToPatch.C
    patchToPatch/intersection/intersectionPatchToPatch.C
    patchToPatch/rays/raysPatchToPatch.C

    patchToPatch/patchToPatchStabilisation/patchToPatchStabilisation.C

    patchToPatchFieldMapper/patchToPatchFieldMapper.C
    patchToPatchFieldMapper/patchToPatchLeftOverFieldMapper.C
    patchToPatchFieldMapper/patchToPatchNormalisedFieldMapper.C

    cellsToCells/cellsToCells/cellsToCells.C
    cellsToCells/cellsToCells/cellsToCellsParallelOps.C
    cellsToCells/matching/matchingCellsToCells.C
    cellsToCells/nearest/nearestCellsToCells.C
    cellsToCells/intersection/intersectionCellsToCells.C

    cellsToCells/cellsToCellsStabilisation/cellsToCellsStabilisation.C

    meshToMesh/meshToMesh.C

    mappedPatches/mappedPatchBaseBase/mappedPatchBaseBase.C

    mappedPatches/mappedPatchBase/mappedPatchBase.C
    mappedPatches/mappedPolyPatch/mappedPolyPatch.C
    mappedPatches/mappedPolyPatch/mappedWallPolyPatch.C
    mappedPatches/mappedPointPatch/mappedPointPatch.C
    mappedPatches/mappedPointPatch/mappedWallPointPatch.C

    mappedPatches/mappedExtrudedPatchBase/mappedExtrudedPatchBase.C
    mappedPatches/mappedExtrudedPolyPatch/mappedExtrudedWallPolyPatch.C
    mappedPatches/mappedExtrudedPointPatch/mappedExtrudedWallPointPatch.C

    mappedPatches/mappedInternalPatchBase/mappedInternalPatchBase.C
    mappedPatches/mappedInternalPolyPatch/mappedInternalPolyPatch.C
    mappedPatches/mappedInternalPointPatch/mappedInternalPointPatch.C

    mappedPatches/nonConformalMappedPatchBase/nonConformalMappedPatchBase.C
    mappedPatches/nonConformalMappedPolyPatch/nonConformalMappedWallPolyPatch.C
    mappedPatches/nonConformalMappedPointPatch/nonConformalMappedWallPointPatch.C

    meshStructure/meshStructure.C
    meshStructure/topoDistanceData.C
    meshStructure/pointTopoDistanceData.C

    tetOverlapVolume/tetOverlapVolume.C

    triIntersect/triIntersect.C
    triIntersect/triIntersectLocationIO.C

    cutPoly/cellEdgeAddressing.C
    cutPoly/cutPoly.C
    cutPoly/cutPolyIsoSurface.C

    nonConformal/polyPatches/nonConformal/nonConformalPolyPatch.C
    nonConformal/polyPatches/nonConformalCoupled/nonConformalCoupledPolyPatch.C
    nonConformal/polyPatches/nonConformalCyclic/nonConformalCyclicPolyPatch.C
    nonConformal/polyPatches/nonConformalError/nonConformalErrorPolyPatch.C
    nonConformal/polyPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPolyPatch.C

    nonConformal/pointPatches/nonConformalCyclic/nonConformalCyclicPointPatch.C
    nonConformal/pointPatches/nonConformalError/nonConformalErrorPointPatch.C
    nonConformal/pointPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPointPatch.C

    nonConformal/pointPatchFields/nonConformalCyclic/nonConformalCyclicPointPatchFields.C
    nonConformal/pointPatchFields/nonConformalError/nonConformalErrorPointPatchFields.C
    nonConformal/pointPatchFields/nonConformalProcessorCyclic/nonConformalProcessorCyclicPointPatchFields.C

    nonConformal/boundary/nonConformalBoundary.C

    pointDist/pointEdgeDist.C
    pointDist/pointDist.C

    patchIntersection/star/star.C
    patchIntersection/patchIntersection.C
    patchIntersection/primitivePatchIntersection.C
    patchIntersection/polyPatchIntersection.C

    patchFaceOrientation/patchFaceOrientation.C

    zoneGenerators/remove/remove.C
    zoneGenerators/clear/clear.C
    zoneGenerators/write/write.C
    zoneGenerators/print/print.C

    zoneGenerators/union/union.C
    zoneGenerators/difference/difference.C
    zoneGenerators/intersection/intersection.C
    zoneGenerators/invert/invert.C

    zoneGenerators/volume/volume/volume.C
    zoneGenerators/volume/box/box.C
    zoneGenerators/volume/sphere/sphere.C
    zoneGenerators/volume/hemisphere/hemisphere.C
    zoneGenerators/volume/cylinder/cylinder.C
    zoneGenerators/volume/annulus/annulus.C
    zoneGenerators/volume/truncatedCone/truncatedCone.C
    zoneGenerators/volume/insideSurface/insideSurface.C

    zoneGenerators/point/point/point_zoneGenerator.C

    zoneGenerators/cell/containsPoints/containsPoints.C
    zoneGenerators/cell/generatedCellZone/generatedCellZone.C
    zoneGenerators/cell/patchCells/patchCells_zoneGenerator.C
    zoneGenerators/cell/nearPatchCells/nearPatchCells_zoneGenerator.C

    zoneGenerators/face/face/face_zoneGenerator.C
    zoneGenerators/face/patch/patch_zoneGenerator.C
    zoneGenerators/face/normal/normal_zoneGenerator.C
    zoneGenerators/face/orient/orient_zoneGenerator.C
    zoneGenerators/face/flip/flip_zoneGenerator.C
    zoneGenerators/face/plane/plane_zoneGenerator.C
    zoneGenerators/face/surface/surface_zoneGenerator.C
    zoneGenerators/face/generatedFaceZone/generatedFaceZone.C

    zoneGenerators/periodic/periodic.C

    zoneGenerators/set/set.C
)

target_include_directories(meshTools PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(meshTools PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
    triSurface
    fileFormats
    Pstream_best
)


target_compile_definitions(meshTools PRIVATE
    LIB_NAME=libmeshTools.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};meshTools" PARENT_SCOPE)


