lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(mixtureViscosityModels)
add_subdirectory(relativeVelocityModels)
add_subdirectory(packingDispersionModels)

set(target_name incompressibleDriftFluxSolver)

add_library(${target_name} SHARED
    incompressibleDriftFluxMixture/incompressibleDriftFluxMixture.C
    alphaSuSp.C
    incompressibleDriftFlux.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${FOAM_MODULES}/VoFSolver/lnInclude
    #    ${FOAM_MODULES}/basicFluidSolver/lnInclude

    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC
    twoPhaseSolver

    incompressibleDriftFluxViscosityModel
    incompressibleDriftFluxRelativeVelocityModels
    incompressibleDriftFluxPackingDispersionModels
    #    physicalProperties
    #    interfaceCompression
    incompressibleTwoPhases
    twoPhaseProperties
    #    interfaceProperties
    momentumTransportModels
    compressibleMomentumTransportModels
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


