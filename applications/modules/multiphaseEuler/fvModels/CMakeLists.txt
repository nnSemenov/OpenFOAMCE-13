lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name multiphaseEulerFvModels)

add_library(${target_name} SHARED
    # Turbulence
    phaseTurbulenceStabilisation/phaseTurbulenceStabilisation.C
    interfaceTurbulenceDamping/interfaceTurbulenceDamping.C

    # Interfacial phase change
    heatTransferLimitedPhaseChange/heatTransferLimitedPhaseChange.C
    massDiffusionLimitedPhaseChange/massDiffusionLimitedPhaseChange.C
    reactionDrivenPhaseChange/reactionDrivenPhaseChange.C

    # Wall phase change
    wallPhaseChange/wallPhaseChange.C
    wallPhaseChange/alphatPhaseChangeWallFunctionFvPatchScalarField.C

    # Nucleation phase change
    nucleation/nucleation.C
    homogeneousNucleation/homogeneousNucleation.C
    homogeneousNucleation/homogeneousCondensation.C
    homogeneousNucleation/homogeneousLiquidPhaseSeparation.C
    wallBoiling/wallBoiling.C
    wallBoiling/wallBoilingPhaseChangeRateFvPatchScalarField.C
    phaseSurfaceBoiling/phaseSurfaceBoiling.C

    # Surface condensation
    wallCondensation/wallCondensation.C
    wallCondensation/wallCondensationPhaseChangeRateFvPatchScalarField.C
    phaseSurfaceCondensation/phaseSurfaceCondensation.C

    # Cavitation
    multiphaseEulerCavitation/multiphaseEulerCavitation.C

    # Source conditions
    derivedFvFieldSources/nucleationGroupFraction/nucleationGroupFractionFvScalarFieldSource.C
    derivedFvFieldSources/nucleationGroupSurfaceAreaVolumeRatio/nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource.C
    derivedFvFieldSources/nucleationInterfacialCurvature/nucleationInterfacialCurvatureFvScalarFieldSource.C

    # Wall boiling sub models
    wallBoiling/partitioningModels/partitioningModel/partitioningModel.C
    wallBoiling/partitioningModels/partitioningModel/partitioningModelNew.C
    wallBoiling/partitioningModels/phaseFraction/phaseFraction.C
    wallBoiling/partitioningModels/Lavieville/Lavieville.C
    wallBoiling/partitioningModels/cosine/cosine.C
    wallBoiling/partitioningModels/linear/linear.C

    wallBoiling/nucleationSiteModels/nucleationSiteModel/nucleationSiteModel.C
    wallBoiling/nucleationSiteModels/nucleationSiteModel/nucleationSiteModelNew.C
    wallBoiling/nucleationSiteModels/LemmertChawla/LemmertChawla.C
    wallBoiling/nucleationSiteModels/KocamustafaogullariIshiiNucleationSite/KocamustafaogullariIshiiNucleationSite.C

    wallBoiling/departureDiameterModels/departureDiameterModel/departureDiameterModel.C
    wallBoiling/departureDiameterModels/departureDiameterModel/departureDiameterModelNew.C
    wallBoiling/departureDiameterModels/TolubinskiKostanchuk/TolubinskiKostanchuk.C
    wallBoiling/departureDiameterModels/KocamustafaogullariIshiiDepartureDiameter/KocamustafaogullariIshiiDepartureDiameter.C

    wallBoiling/departureFrequencyModels/departureFrequencyModel/departureFrequencyModel.C
    wallBoiling/departureFrequencyModels/departureFrequencyModel/departureFrequencyModelNew.C
    wallBoiling/departureFrequencyModels/Cole/Cole.C
    wallBoiling/departureFrequencyModels/KocamustafaogullariIshiiDepartureFrequency/KocamustafaogullariIshiiDepartureFrequency.C

    # Sintering models
    KochFriedlanderSintering/KochFriedlanderSintering.C

    # Capillary models
    capillary/capillaryForce.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${FOAM_MODULES}/basicFluidSolver/lnInclude
    #    ${LIB_SRC}/meshTools/lnInclude
    #    ${LIB_SRC}/finiteVolume/lnInclude
    #    ${LIB_SRC}/sampling/lnInclude
    #    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude
    ${LIB_SRC}/thermophysicalModels/saturationModels/lnInclude
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    ${LIB_SRC}/twoPhaseModels/compressibleCavitation/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/coupledThermophysicalTransportModels/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)

target_link_libraries(${target_name} PUBLIC
    multiphaseEulerSolver

    fvModels
    phaseSystem
    populationBalance
    compressibleCavitationModels
    phaseFluidThermophysicalTransportModels
    phaseFluidMulticomponentThermophysicalTransportModels
    coupledThermophysicalTransportModels
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

target_compile_features(${target_name} PRIVATE
    cxx_std_23
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


