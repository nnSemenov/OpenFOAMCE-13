lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name sixDoFRigidBodyMeshMotion)

#wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)

add_library(${target_name} SHARED
#    ${sources_from_wmake}
    sixDoFRigidBodyMeshMotion.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
#    ${LIB_SRC}/rigidBodyMotion/sixDoFRigidBodyMotion/lnInclude
    ${LIB_SRC}/finiteVolume/lnInclude 
#    ${LIB_SRC}/meshTools/lnInclude
#    ${LIB_SRC}/functionObjects/forces/lnInclude
#    ${LIB_SRC}/fileFormats/lnInclude
    ${LIB_SRC}/motionSolvers/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    sixDoFRigidBodyMotion
    forces
    meshTools
    fileFormats
    motionSolver_fvMeshMover
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


