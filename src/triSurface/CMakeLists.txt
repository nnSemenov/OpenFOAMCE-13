lnInclude(${CMAKE_CURRENT_SOURCE_DIR})


set(triSurfaceTools       triSurface/triSurfaceTools)
set(geometricSurfacePatch triSurface/geometricSurfacePatch)
set(interfaces            triSurface/interfaces)

add_library(triSurface SHARED

    meshTriangulation/meshTriangulation.C

    triSurface/triSurface.C
    triSurface/triSurfaceAddressing.C
    triSurface/stitchTriangles.C

    ${interfaces}/STL/writeSTL.C
    ${interfaces}/STL/readSTL.C
    #${interfaces}/STL/readSTLASCII.L
    ${interfaces}/STL/readSTLBINARY.C
    ${interfaces}/GTS/writeGTS.C
    ${interfaces}/GTS/readGTS.C
    ${interfaces}/OBJ/readOBJ.C
    ${interfaces}/OBJ/writeOBJ.C
    ${interfaces}/SMESH/writeSMESH.C
    ${interfaces}/OFF/readOFF.C
    ${interfaces}/OFF/writeOFF.C
    ${interfaces}/TRI/writeTRI.C
    ${interfaces}/TRI/readTRI.C
    ${interfaces}/AC3D/readAC.C
    ${interfaces}/AC3D/writeAC.C
    ${interfaces}/VTK/readVTK.C
    ${interfaces}/VTK/writeVTK.C
    ${interfaces}/NAS/readNAS.C

    triSurface/geometricSurfacePatch/geometricSurfacePatch.C
    triSurface/surfacePatch/surfacePatch.C
    triSurface/surfacePatch/surfacePatchIOList.C

    tools/labelledTri/sortLabelledTri.C
)

flex_target(triSurfaceFlex
    ${interfaces}/STL/readSTLASCII.L
    realSTLASCII.C
    COMPILE_FLAGS "-+ -f"
)
target_sources(triSurface PRIVATE
    ${FLEX_triSurfaceFlex_OUTPUTS}
)

target_include_directories(triSurface PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(triSurface PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
    fileFormats
    surfMesh
)

target_compile_definitions(triSurface PRIVATE
    LIB_NAME=libtriSurface.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};triSurface" PARENT_SCOPE)


