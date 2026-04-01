lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(surfMesh SHARED
    surfZone/surfZone/surfZone.C
    surfZone/surfZone/surfZoneIOList.C
    surfZone/surfZoneIdentifier/surfZoneIdentifier.C

    MeshedSurfaceAllocator/MeshedSurfaceIOAllocator.C

    MeshedSurface/MeshedSurfaceCore.C
    MeshedSurface/MeshedSurfaces.C
    UnsortedMeshedSurface/UnsortedMeshedSurfaces.C

    MeshedSurfaceProxy/MeshedSurfaceProxyCore.C

    surfaceRegistry/surfaceRegistry.C
    surfMesh/surfMesh.C
    surfMesh/surfMeshClear.C
    surfMesh/surfMeshIO.C

    surfaceFormats/surfaceFormatsCore.C

    surfaceFormats/ac3d/AC3DsurfaceFormatCore.C
    surfaceFormats/ac3d/AC3DsurfaceFormatRunTime.C
    surfaceFormats/gts/GTSsurfaceFormatRunTime.C
    surfaceFormats/nas/NASsurfaceFormatRunTime.C
    surfaceFormats/obj/OBJsurfaceFormatRunTime.C
    surfaceFormats/off/OFFsurfaceFormatRunTime.C
    surfaceFormats/ofs/OFSsurfaceFormatCore.C
    surfaceFormats/ofs/OFSsurfaceFormatRunTime.C
    surfaceFormats/smesh/SMESHsurfaceFormatRunTime.C
    surfaceFormats/starcd/STARCDsurfaceFormatCore.C
    surfaceFormats/starcd/STARCDsurfaceFormatRunTime.C
    surfaceFormats/stl/STLsurfaceFormatCore.C
    surfaceFormats/stl/STLsurfaceFormatRunTime.C
    #surfaceFormats/stl/STLsurfaceFormatASCII.L
    surfaceFormats/tri/TRIsurfaceFormatCore.C
    surfaceFormats/tri/TRIsurfaceFormatRunTime.C
    surfaceFormats/vtk/VTKsurfaceFormatCore.C
    surfaceFormats/vtk/VTKsurfaceFormatRunTime.C
    surfaceFormats/wrl/WRLsurfaceFormatCore.C
    surfaceFormats/wrl/WRLsurfaceFormatRunTime.C
    surfaceFormats/x3d/X3DsurfaceFormatCore.C
    surfaceFormats/x3d/X3DsurfaceFormatRunTime.C
)

flex_target(surfMeshFlex
    surfaceFormats/stl/STLsurfaceFormatASCII.L
    STLsurfaceFormatASCII.C
    COMPILE_FLAGS "-+ -f"
)
target_sources(surfMesh PRIVATE
    ${FLEX_surfMeshFlex_OUTPUTS}
)

target_include_directories(surfMesh PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_compile_definitions(surfMesh PRIVATE
    LIB_NAME=libsurfMesh.so
)

target_link_libraries(surfMesh PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
    fileFormats
#    FLEX::fl
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};surfMesh" PARENT_SCOPE)


