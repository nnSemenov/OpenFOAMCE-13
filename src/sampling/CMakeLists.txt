lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(sampling SHARED
    probes/probes.C
    probes/patchProbes.C
    probes/probesGrouping.C

    coordSet/coordSet.C

    sampledSet/sampledSet/sampledSet.C
    sampledSet/sampledSet/sampledSetParticle.C
    sampledSet/sampledSet/sampledSetCloud.C
    sampledSet/sampledSets/sampledSets.C

    sampledSet/arcUniform/arcUniform.C
    sampledSet/boxUniform/boxUniform.C
    sampledSet/circleRandom/circleRandom.C
    sampledSet/lineCell/lineCell.C
    sampledSet/lineCellFace/lineCellFace.C
    sampledSet/lineFace/lineFace.C
    sampledSet/lineUniform/lineUniform.C
    sampledSet/boundaryPoints/boundaryPoints.C
    sampledSet/boundaryRandom/boundaryRandom.C
    sampledSet/points/points.C
    sampledSet/sphereRandom/sphereRandom.C
    sampledSet/triSurface/triSurfaceSampledSet.C
    sampledSet/cellSet/cellSetSampledSet.C
    sampledSet/faceSet/faceSetSampledSet.C
)

set(setWriters sampledSet/writers)
target_sources(sampling PRIVATE
    ${setWriters}/setWriter.C

    ${setWriters}/ensight/file/ensightFile.C
    ${setWriters}/ensight/file/ensightGeoFile.C
    ${setWriters}/ensight/part/ensightPart.C
    ${setWriters}/ensight/part/ensightPartIO.C
    ${setWriters}/ensight/part/ensightPartCells.C
    ${setWriters}/ensight/part/ensightPartFaces.C
    ${setWriters}/ensight/part/ensightParts.C
    ${setWriters}/ensight/ensightSetWriter.C
    ${setWriters}/raw/rawSetWriter.C
    ${setWriters}/vtk/vtkSetWriter.C
    ${setWriters}/csv/csvSetWriter.C
    ${setWriters}/gnuplot/gnuplotSetWriter.C
    ${setWriters}/none/noSetWriter.C
)

target_sources(sampling PRIVATE
    sampledSurface/sampledPatch/sampledPatch.C
    sampledSurface/sampledPatchInternalField/sampledPatchInternalField.C
    sampledSurface/sampledIsoSurfaceSurface/sampledIsoSurfaceSurface.C
    sampledSurface/sampledIsoSurface/sampledIsoSurface.C
    sampledSurface/sampledDistanceSurface/sampledDistanceSurface.C
    sampledSurface/sampledCutPlane/sampledCutPlane.C
    sampledSurface/sampledSurface/sampledSurface.C
    sampledSurface/sampledSurfaces/sampledSurfaces.C
    sampledSurface/sampledTriSurface/sampledTriSurface.C
    sampledSurface/sampledThresholdCellFaces/thresholdCellFaces.C
    sampledSurface/sampledThresholdCellFaces/sampledThresholdCellFaces.C
)

set(surfWriters sampledSurface/writers)
target_sources(sampling PRIVATE
    ${surfWriters}/surfaceWriter.C
    ${surfWriters}/none/noSurfaceWriter.C
    ${surfWriters}/ensight/ensightSurfaceWriter.C
    ${surfWriters}/ensight/ensightPTraits.C
    ${surfWriters}/foam/foamSurfaceWriter.C
    ${surfWriters}/proxy/proxySurfaceWriter.C
    ${surfWriters}/raw/rawSurfaceWriter.C
    ${surfWriters}/vtk/vtkSurfaceWriter.C
)

target_sources(sampling PRIVATE
    cutPlot/patchCutPlot.C
    cutPlot/cellCutPlot.C
)

target_include_directories(sampling PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(sampling PRIVATE
    ${LIB_SRC}/tracking/lnInclude
)

target_link_libraries(sampling PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    finiteVolume 
    meshTools 
    surfMesh 
    fileFormats 
    triSurface 
    lagrangian
)

target_compile_definitions(sampling PRIVATE
    LIB_NAME=libsampling.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};sampling" PARENT_SCOPE)


