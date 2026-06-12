lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name rigidBodyPropellerDisk)

add_library(${target_name} SHARED
    rigidBodyPropellerDisk.C
    propellerDiskForce/propellerDiskForce.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/functionObjects/forces/lnInclude
    ${LIB_SRC}/finiteVolume/lnInclude
    ${LIB_SRC}/meshTools/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PRIVATE
    propellerDisk
    pointMeshMovers
    rigidBodyDynamics
    rigidBodyMeshMotion
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


