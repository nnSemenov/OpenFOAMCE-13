lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

find_package(VTK
    COMPONENTS CommonCore
    REQUIRED
)

find_package(ParaView REQUIRED)

set(target_name vtkPVFoam)

add_library(${target_name} SHARED
    vtkPVFoam.C
    vtkPVFoamFields.C
    vtkPVFoamMesh.C
    vtkPVFoamMeshLagrangian.C
    vtkPVFoamMeshSet.C
    vtkPVFoamMeshVolume.C
    vtkPVFoamMeshZone.C
    vtkPVFoamUpdateInfo.C
    vtkPVFoamUtils.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ../PVFoamReader/vtk
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    finiteVolume
    meshTools
    parallel
    polyTopoChange
    genericFvFields
    genericLagrangianFields
    lagrangian
    Lagrangian

    VTK::CommonCore
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


