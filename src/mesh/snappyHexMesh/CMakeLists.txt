lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name snappyHexMesh)

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)

add_library(${target_name} SHARED
    ${sources_from_wmake}
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/parallel/decompose/decompositionMethods/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC
    polyTopoChange
    pointMeshMovers 
    finiteVolume 
    lagrangian 
    meshTools 
    meshCheck 
    fileFormats 
    surfMesh 
    triSurface 
    distributed
    tracking
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


#polyTopoChange
#motionSolvers
#finiteVolume
#lagrangian
#meshTools
#meshCheck
#fileFormats
#surfMesh
#triSurface
#distributed