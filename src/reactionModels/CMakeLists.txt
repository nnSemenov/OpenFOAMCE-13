lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name reactionModels)

add_library(${target_name} SHARED
    reactionModel/reactionModel.C
    reactionModel/reactionModelNew.C

    singleStepReaction/singleStepReaction.C

    diffusion/diffusion.C

    infinitelyFastChemistry/infinitelyFastChemistry.C

    PaSR/PaSR.C

    EDC/EDC.C

    laminar/laminar.C

    FSD/reactionRateFlameAreaModels/consumptionSpeed/consumptionSpeed.C
    FSD/reactionRateFlameAreaModels/reactionRateFlameArea/reactionRateFlameArea.C
    FSD/reactionRateFlameAreaModels/reactionRateFlameArea/reactionRateFlameAreaNew.C
    FSD/reactionRateFlameAreaModels/relaxation/relaxation.C
    FSD/FSD.C

    zoneReaction/zoneReaction.C

    noReaction/noReaction.C

    functionObjects/adjustTimeStepToReaction/adjustTimeStepToReaction.C
    functionObjects/Qdot/Qdot.C

    radiationModels/absorptionEmissionModels/greyMeanReaction/greyMeanReaction.C
    radiationModels/absorptionEmissionModels/wideBandReaction/wideBandReaction.C
    radiationModels/sootModels/mixtureFraction/mixtureFraction.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    momentumTransportModels 
    compressibleMomentumTransportModels 
    chemistryModel 
    radiationModels 
    finiteVolume 
    meshTools
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


