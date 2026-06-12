lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

find_package(VTK
    COMPONENTS CommonCore
    REQUIRED)

set(target_name vtkPVblockMesh)

add_library(${target_name} SHARED
    vtkPVblockMesh.C
    vtkPVblockMeshConvert.C
    vtkPVblockMeshUtils.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ../PVblockMeshReader/vtk
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    meshTools
    fileFormats
    blockMesh

    VTK::CommonCore
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


