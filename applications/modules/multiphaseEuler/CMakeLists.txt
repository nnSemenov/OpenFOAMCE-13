lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(phaseSystem)
add_subdirectory(populationBalance)
add_subdirectory(momentumTransportModels)
add_subdirectory(thermophysicalTransportModels)
add_subdirectory(reactions)

set(target_name multiphaseEulerSolver)

add_library(${target_name} SHARED
    setRDeltaT.C
    moveMesh.C
    momentumPredictor.C
    thermophysicalPredictor.C
    compressibilityEqns.C
    facePressureCorrector.C
    cellPressureCorrector.C
    pressureCorrector.C
    multiphaseEuler.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/meshTools/lnInclude
    ${LIB_SRC}/finiteVolume/lnInclude
    ${LIB_SRC}/sampling/lnInclude
    #${LIB_SRC}/physicalProperties/lnInclude
    #${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    basicFluidSolver

    phaseSystem
    populationBalance
    physicalProperties
    fluidThermophysicalModels
    fluidThermophysicalTransportModel
    multiphaseEulerMomentumTransportModels
    multiphaseEulerThermophysicalTransportModels
    multiphaseEulerReactions
    fvModels
    fvConstraints
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

add_subdirectory(fvModels)
add_subdirectory(functionObjects)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)
