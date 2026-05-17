lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name incompressibleVoFSolver)

add_library(${target_name} SHARED
    incompressibleTwoPhaseVoFMixture/incompressibleTwoPhaseVoFMixture.C
    incompressibleInterPhaseTransportModel/incompressibleInterPhaseTransportModel.C
    alphaSuSp.C
    incompressibleVoF.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${FOAM_MODULES}/twoPhaseSolver/lnInclude
    ${LIB_SRC}/twoPhaseModels/interfaceCompression/lnInclude
    ${LIB_SRC}/twoPhaseModels/interfaceProperties/lnInclude
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    ${LIB_SRC}/twoPhaseModels/VoF
    ${LIB_SRC}/twoPhaseModels/incompressibleTwoPhases/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseIncompressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/incompressible/lnInclude
    #    ${FOAM_MODULES}/VoFSolver/lnInclude
    #    ${FOAM_MODULES}/basicFluidSolver/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)

target_link_libraries(${target_name} PUBLIC
    twoPhaseVoFSolver

    #    physicalProperties
    #    interfaceCompression
    incompressibleTwoPhases
    twoPhaseProperties
    #    interfaceProperties
    #    momentumTransportModels
    incompressibleMomentumTransportModels
    phaseIncompressibleMomentumTransportModels
    #    finiteVolume
    #    meshTools
    #    fvModels
    #    fvConstraints
    #    sampling
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}")

add_subdirectory(fvModels/VoFCavitation)
add_subdirectory(fvModels/VoFTurbulenceDamping)

set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)