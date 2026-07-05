lnInclude(${CMAKE_CURRENT_SOURCE_DIR}
    EXCLUDE_DIRS
    fvMeshTopoChangers
    meshToMes
)

set(target_name fvMeshTopoChangers)

add_library(${target_name} SHARED
    refiner/refiner_fvMeshTopoChanger.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    triSurface 
    meshTools 
    polyTopoChange 
    finiteVolume 
    fvMeshStitchers
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


add_subdirectory(meshToMesh)
set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)