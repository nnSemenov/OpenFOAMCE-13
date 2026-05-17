lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name twoPhaseVoFSolver)

add_library(${target_name} SHARED
    setInterfaceRDeltaT.C
    twoPhaseVoFSolver.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${FOAM_MODULES}/basicFluidSolver/lnInclude
    #    ${FOAM_MODULES}/VoFSolver/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC
    twoPhaseSolver

    #    physicalProperties
    #    interfaceCompression
    #    interfaceProperties
    #    twoPhaseMixture
    #    finiteVolume
    #    meshTools
    #    fvModels
    #    fvConstraints
    #    sampling
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


