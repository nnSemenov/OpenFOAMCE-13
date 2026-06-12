lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name populationBalance)

add_library(${target_name} SHARED
    diameterModels/populationBalance/populationBalance.C

    populationBalanceModel/populationBalanceModel.C

    oneDimensionalDiscretisations/oneDimensionalDiscretisation/oneDimensionalDiscretisation.C
    oneDimensionalDiscretisations/uniform/uniform_oneDimensionalDiscretisation.C
    oneDimensionalDiscretisations/exponential/exponential_oneDimensionalDiscretisation.C
    oneDimensionalDiscretisations/manual/manual_oneDimensionalDiscretisation.C

    shapeModels/shapeModel/shapeModel.C
    shapeModels/spherical/spherical.C
    shapeModels/fractal/fractal.C

    coalescenceModels/coalescenceModel/coalescenceModel.C
    coalescenceModels/noCoalescence/noCoalescence.C
    coalescenceModels/AdachiStuartFokkink/AdachiStuartFokkink.C
    coalescenceModels/ballisticCollisions/ballisticCollisions.C
    coalescenceModels/BrownianCollisions/BrownianCollisions.C
    coalescenceModels/constantCoalescence/constantCoalescence.C
    coalescenceModels/CoulaloglouTavlarides/CoulaloglouTavlarides.C
    coalescenceModels/DahnekeInterpolation/DahnekeInterpolation.C
    coalescenceModels/hydrodynamic/hydrodynamic.C
    coalescenceModels/LehrMilliesMewesCoalescence/LehrMilliesMewesCoalescence.C
    coalescenceModels/LiaoCoalescence/LiaoCoalescence.C
    coalescenceModels/Luo/Luo.C
    coalescenceModels/PrinceBlanch/PrinceBlanch.C
    coalescenceModels/turbulentShear/turbulentShear.C

    breakupModels/breakupModel/breakupModel.C
    breakupModels/noBreakup/noBreakup.C
    breakupModels/daughterSizeDistribution/daughterSizeDistribution.C
    breakupModels/exponential/exponential.C
    breakupModels/Kusters/Kusters.C
    breakupModels/Laakkonen/Laakkonen.C
    breakupModels/powerLaw/powerLaw.C
    breakupModels/LehrMilliesMewes/LehrMilliesMewes.C
    breakupModels/Liao/LiaoBase.C
    breakupModels/Liao/Liao.C
    breakupModels/LuoSvendsen/LuoSvendsen.C
    breakupModels/powerLawUniformBinary/powerLawUniformBinary.C

    daughterSizeDistributionModels/daughterSizeDistributionModel/daughterSizeDistributionModel.C
    daughterSizeDistributionModels/LaakkonenDaughterSizeDistribution/LaakkonenDaughterSizeDistribution.C
    daughterSizeDistributionModels/uniformBinary/uniformBinary.C

    derivedFvFields/groupProperty/groupPropertyFvScalarField.C

    derivedFvPatchFields/singleGroupFraction/singleGroupFractionFvPatchScalarField.C
    derivedFvPatchFields/distributionGroupFraction/distributionGroupFractionFvPatchScalarField.C

    derivedFvFieldSources/growth/growthFvScalarFieldSource.C

    derivedFvFieldSources/singleGroupFraction/singleGroupFractionFvScalarFieldSource.C
    derivedFvFieldSources/distributionGroupFraction/distributionGroupFractionFvScalarFieldSource.C
    derivedFvFieldSources/growthGroupFraction/growthGroupFractionFvScalarFieldSource.C
    derivedFvFieldSources/interfacialGrowthGroupFraction/interfacialGrowthGroupFractionFvScalarFieldSource.C

    derivedFvFieldSources/growthSecondaryProperty/growthSecondaryPropertyFvScalarFieldSource.C

    derivedFvFieldSources/uniformFixedValueGroupSurfaceAreaVolumeRatio/uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource.C
    derivedFvFieldSources/interfacialGrowthGroupSurfaceAreaVolumeRatio/interfacialGrowthGroupSurfaceAreaVolumeRatioFvScalarFieldSource.C
    derivedFvFieldSources/hardSphereGroupSurfaceAreaVolumeRatio/hardSphereGroupSurfaceAreaVolumeRatioFvScalarFieldSource.C
    derivedFvFieldSources/ParkRogakGroupSurfaceAreaVolumeRatio/ParkRogakGroupSurfaceAreaVolumeRatioFvScalarFieldSource.C

    populationBalanceSystem/populationBalanceSystem.C
)


target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${LIB_SRC}/meshTools/lnInclude
    #    ${LIB_SRC}/finiteVolume/lnInclude
    #    ${LIB_SRC}/sampling/lnInclude
    #    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude
    #    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    #    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC
    phaseSystem
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


