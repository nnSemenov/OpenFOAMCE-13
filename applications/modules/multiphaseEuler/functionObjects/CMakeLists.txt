lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name multiphaseEulerFunctionObjects)

add_library(${target_name} SHARED
    populationBalanceMoments/populationBalanceMoments.C
    populationBalanceSizeDistribution/populationBalanceSizeDistribution.C
    populationBalanceSetSizeDistribution/populationBalanceSetSizeDistribution.C

    phaseForces/phaseForces.C
    phaseMap/phaseMap.C

    wallBoilingProperty/wallBoilingProperty.C

    adjustTimeStepToNucleation/adjustTimeStepToNucleation.C

    capillary/capillary.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${LIB_SRC}/meshTools/lnInclude
    #    ${LIB_SRC}/finiteVolume/lnInclude
    #    ${LIB_SRC}/sampling/lnInclude
    ${LIB_SRC}/functionObjects/field/lnInclude
    #    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
    ${LIB_SRC}/fvModels/general/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC
    phaseSystem
    populationBalance
    multiphaseEulerFvModels
    fileFormats
    fieldFunctionObjects
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


