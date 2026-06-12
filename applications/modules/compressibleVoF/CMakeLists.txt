lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name compressibleVoFSolver)

add_library(${target_name} SHARED
    compressibleTwoPhaseVoFMixture/compressibleTwoPhaseVoFMixture.C
    compressibleInterPhaseTransportModel/compressibleInterPhaseTransportModel.C
    compressibleInterPhaseThermophysicalTransportModel/compressibleInterPhaseThermophysicalTransportModel.C

    alphaSuSp.C
    momentumPredictor.C
    thermophysicalPredictor.C
    pressureCorrector.C
    compressibleVoF.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/twoPhaseModels/VoF
    ${LIB_SRC}/twoPhaseModels/interfaceCompression/lnInclude
    ${LIB_SRC}/twoPhaseModels/interfaceProperties/lnInclude
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude

    #    ${FOAM_MODULES}/twoPhaseSolver/lnInclude
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
    fluidThermophysicalModels
    #    interfaceCompression
    compressibleInterfaceProperties
    compressibleTwoPhases
    momentumTransportModels
    compressibleMomentumTransportModels
    phaseCompressibleMomentumTransportModels
    thermophysicalTransportModel
    coupledThermophysicalTransportModels
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

add_subdirectory(fvModels/VoFTurbulenceDamping)
add_subdirectory(fvModels/VoFClouds)
add_subdirectory(fvModels/VoFCavitation)
add_subdirectory(fvModels/VoFSolidificationMelting)

set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)