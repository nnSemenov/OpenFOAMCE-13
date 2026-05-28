
find_package(ParaView REQUIRED)

add_subdirectory(vtkPVblockMesh)
add_subdirectory(vtkPVFoam)


set(BUILD_SHARED_LIBS ON)

#add_subdirectory(PVFoamReader)

PARAVIEW_PLUGIN_SCAN(
    PLUGIN_FILES
        "PVblockMeshReader/paraview.plugin"
        "PVFoamReader/paraview.plugin"
    PROVIDES_PLUGINS plugins
    ENABLE_BY_DEFAULT ON
)

PARAVIEW_PLUGIN_BUILD(PLUGINS ${plugins})

set(FOAM_regular_executables ${FOAM_regular_executables} PARENT_SCOPE)
set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)
set(FOAM_special_libraries "${FOAM_special_libraries}" PARENT_SCOPE)