lnInclude(${CMAKE_CURRENT_SOURCE_DIR} HEADER_LIST target_headers)

set(target_name phaseSystem)

add_library(${target_name} SHARED
    ${target_headers}

    phaseModel/phaseModel.C
    phaseModel/phaseModelNew.C
    phaseModels/phaseModels.C

    phaseInterface/phaseInterfaceKey/phaseInterfaceKey.C
    phaseInterface/phaseInterface/phaseInterface.C
    phaseInterface/dispersedPhaseInterface/dispersedPhaseInterface.C
    phaseInterface/segregatedPhaseInterface/segregatedPhaseInterface.C
    phaseInterface/displacedPhaseInterface/displacedPhaseInterface.C
    phaseInterface/dispersedDisplacedPhaseInterface/dispersedDisplacedPhaseInterface.C
    phaseInterface/segregatedDisplacedPhaseInterface/segregatedDisplacedPhaseInterface.C
    phaseInterface/sidedPhaseInterface/sidedPhaseInterface.C
    phaseInterface/dispersedSidedPhaseInterface/dispersedSidedPhaseInterface.C
    phaseInterface/segregatedSidedPhaseInterface/segregatedSidedPhaseInterface.C
    phaseInterface/displacedSidedPhaseInterface/displacedSidedPhaseInterface.C
    phaseInterface/dispersedDisplacedSidedPhaseInterface/dispersedDisplacedSidedPhaseInterface.C
    phaseInterface/segregatedDisplacedSidedPhaseInterface/segregatedDisplacedSidedPhaseInterface.C

    phaseSystem/phaseSystem.C
    phaseSystem/phaseSystemSolve.C

    diameterModels/diameterModel/diameterModel.C
    diameterModels/diameterModel/diameterModelNew.C
    diameterModels/sphericalDiameter/sphericalDiameter.C
    diameterModels/constantDiameter/constantDiameter.C
    diameterModels/isothermalDiameter/isothermalDiameter.C
    diameterModels/linearTsubDiameter/linearTsubDiameter.C
    diameterModels/noDiameter/noDiameter.C
    diameterModels/residualDiameter/residualDiameter.C
    diameterModels/fixedInterfacialAreaDiameter/fixedInterfacialAreaDiameter.C

    diameterModels/IATE/IATE.C
    diameterModels/IATE/IATEsources/IATEsource/IATEsource.C
    diameterModels/IATE/IATEsources/wakeEntrainmentCoalescence/IATEwakeEntrainmentCoalescence.C
    diameterModels/IATE/IATEsources/turbulentBreakUp/IATEturbulentBreakUp.C
    diameterModels/IATE/IATEsources/randomCoalescence/IATErandomCoalescence.C
    diameterModels/IATE/derivedFvFieldSources/interfacialGrowthInterfacialCurvature/interfacialGrowthInterfacialCurvatureFvScalarFieldSource.C

    surfaceTensionCoefficientModels/surfaceTensionCoefficientModel/surfaceTensionCoefficientModel.C
    surfaceTensionCoefficientModels/surfaceTensionCoefficientModel/surfaceTensionCoefficientModelNew.C
    surfaceTensionCoefficientModels/constantSurfaceTensionCoefficient/constantSurfaceTensionCoefficient.C

    interfacialModels/BlendedInterfacialModel/blendingMethods/blendingMethod/blendingMethod.C
    interfacialModels/BlendedInterfacialModel/blendingMethods/blendingMethod/blendingMethodNew.C
    interfacialModels/BlendedInterfacialModel/blendingMethods/linear/linear.C
    interfacialModels/BlendedInterfacialModel/blendingMethods/hyperbolic/hyperbolic.C
    interfacialModels/BlendedInterfacialModel/blendingMethods/continuous/continuous.C
    interfacialModels/BlendedInterfacialModel/blendingMethods/segregatedBlendingMethod/segregatedBlendingMethod.C

    interfacialModels/dragModels/dragModel/dragModel.C
    interfacialModels/dragModels/dragModel/dragModelNew.C
    interfacialModels/dragModels/dispersedDragModel/dispersedDragModel.C
    interfacialModels/dragModels/Beetstra/Beetstra.C
    interfacialModels/dragModels/segregated/segregated.C
    interfacialModels/dragModels/Ergun/Ergun.C
    interfacialModels/dragModels/Gibilaro/Gibilaro.C
    interfacialModels/dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C
    interfacialModels/dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C
    interfacialModels/dragModels/Lain/Lain.C
    interfacialModels/dragModels/SchillerNaumann/SchillerNaumann.C
    interfacialModels/dragModels/SyamlalOBrien/SyamlalOBrien.C
    interfacialModels/dragModels/TomiyamaCorrelated/TomiyamaCorrelated.C
    interfacialModels/dragModels/TomiyamaAnalytic/TomiyamaAnalytic.C
    interfacialModels/dragModels/Tenneti/Tenneti.C
    interfacialModels/dragModels/TomiyamaKataokaZunSakaguchi/TomiyamaKataokaZunSakaguchi.C
    interfacialModels/dragModels/WenYu/WenYu.C
    interfacialModels/dragModels/IshiiZuber/IshiiZuber.C
    interfacialModels/dragModels/AttouFerschneider/AttouFerschneider.C
    interfacialModels/dragModels/aerosolDrag/aerosolDrag.C
    interfacialModels/dragModels/timeScaleFilteredDrag/timeScaleFilteredDrag.C

    interfacialModels/swarmCorrections/swarmCorrection/swarmCorrection.C
    interfacialModels/swarmCorrections/swarmCorrection/swarmCorrectionNew.C
    interfacialModels/swarmCorrections/noSwarm/noSwarm.C
    interfacialModels/swarmCorrections/TomiyamaSwarm/TomiyamaSwarm.C

    interfacialModels/liftModels/liftModel/liftModel.C
    interfacialModels/liftModels/liftModel/liftModelNew.C
    interfacialModels/liftModels/dispersedLiftModel/dispersedLiftModel.C
    interfacialModels/liftModels/noLift/noLift.C
    interfacialModels/liftModels/constantLiftCoefficient/constantLiftCoefficient.C
    interfacialModels/liftModels/Moraga/Moraga.C
    interfacialModels/liftModels/LegendreMagnaudet/LegendreMagnaudet.C
    interfacialModels/liftModels/SaffmanMei/SaffmanMei.C
    interfacialModels/liftModels/TomiyamaLift/TomiyamaLift.C
    interfacialModels/liftModels/wallDampedLift/wallDampedLift.C

    interfacialModels/virtualMassModels/virtualMassModel/virtualMassModel.C
    interfacialModels/virtualMassModels/virtualMassModel/virtualMassModelNew.C
    interfacialModels/virtualMassModels/dispersedVirtualMassModel/dispersedVirtualMassModel.C
    interfacialModels/virtualMassModels/noVirtualMass/noVirtualMass.C
    interfacialModels/virtualMassModels/constantVirtualMassCoefficient/constantVirtualMassCoefficient.C
    interfacialModels/virtualMassModels/Lamb/Lamb.C

    interfacialModels/wallLubricationModels/wallLubricationModel/wallLubricationModel.C
    interfacialModels/wallLubricationModels/wallLubricationModel/wallLubricationModelNew.C
    interfacialModels/wallLubricationModels/dispersedWallLubricationModel/dispersedWallLubricationModel.C
    interfacialModels/wallLubricationModels/noWallLubrication/noWallLubrication.C
    interfacialModels/wallLubricationModels/Antal/Antal.C
    interfacialModels/wallLubricationModels/Frank/Frank.C
    interfacialModels/wallLubricationModels/TomiyamaWallLubrication/TomiyamaWallLubrication.C

    interfacialModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.C
    interfacialModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModelNew.C
    interfacialModels/turbulentDispersionModels/dispersedTurbulentDispersionModel/dispersedTurbulentDispersionModel.C
    interfacialModels/turbulentDispersionModels/noTurbulentDispersion/noTurbulentDispersion.C
    interfacialModels/turbulentDispersionModels/constantTurbulentDispersionCoefficient/constantTurbulentDispersionCoefficient.C
    interfacialModels/turbulentDispersionModels/Burns/Burns.C
    interfacialModels/turbulentDispersionModels/Gosman/Gosman.C
    interfacialModels/turbulentDispersionModels/LopezDeBertodano/LopezDeBertodano.C

    interfacialModels/aspectRatioModels/aspectRatioModel/aspectRatioModel.C
    interfacialModels/aspectRatioModels/aspectRatioModel/aspectRatioModelNew.C
    interfacialModels/aspectRatioModels/constantAspectRatio/constantAspectRatio.C
    interfacialModels/aspectRatioModels/TomiyamaAspectRatio/TomiyamaAspectRatio.C
    interfacialModels/aspectRatioModels/VakhrushevEfremov/VakhrushevEfremov.C
    interfacialModels/aspectRatioModels/Wellek/Wellek.C

    interfacialModels/wallDependentModel/wallDependentModel.C

    interfacialModels/wallDampingModels/wallDampingModel/wallDampingModel.C
    interfacialModels/wallDampingModels/wallDampingModel/wallDampingModelNew.C
    interfacialModels/wallDampingModels/linear/linearWallDamping.C
    interfacialModels/wallDampingModels/cosine/cosineWallDamping.C
    interfacialModels/wallDampingModels/sine/sineWallDamping.C

    interfacialModels/heatTransferModels/heatTransferModel/heatTransferModel.C
    interfacialModels/heatTransferModels/heatTransferModel/heatTransferModelNew.C
    interfacialModels/heatTransferModels/constantNu/constantNuHeatTransfer.C
    interfacialModels/heatTransferModels/RanzMarshall/RanzMarshall.C
    interfacialModels/heatTransferModels/Gunn/Gunn.C
    interfacialModels/heatTransferModels/sphericalHeatTransfer/sphericalHeatTransfer.C
    interfacialModels/heatTransferModels/nonSphericalHeatTransfer/nonSphericalHeatTransfer.C
    interfacialModels/heatTransferModels/timeScaleFilteredHeatTransfer/timeScaleFilteredHeatTransfer.C
    interfacialModels/heatTransferModels/Prandtl/Prandtl.C

    interfacialModels/diffusiveMassTransferModels/diffusiveMassTransferModel/diffusiveMassTransferModel.C
    interfacialModels/diffusiveMassTransferModels/diffusiveMassTransferModel/diffusiveMassTransferModelNew.C
    interfacialModels/diffusiveMassTransferModels/Frossling/Frossling.C
    interfacialModels/diffusiveMassTransferModels/sphericalDiffusiveMassTransfer/sphericalDiffusiveMassTransfer.C

    interfacialModels/interfaceCompositionModels/interfaceCompositionModel/interfaceCompositionModel.C
    interfacialModels/interfaceCompositionModels/interfaceCompositionModel/interfaceCompositionModelNew.C
    interfacialModels/interfaceCompositionModels/Henry/Henry.C
    interfacialModels/interfaceCompositionModels/nonRandomTwoLiquid/nonRandomTwoLiquid.C
    interfacialModels/interfaceCompositionModels/Raoult/Raoult.C
    interfacialModels/interfaceCompositionModels/saturated/saturated.C

    momentumTransferSystem/momentumTransferSystem.C
    heatTransferSystem/heatTransferSystem.C

    oneResistanceHeatTransfer/oneResistanceHeatTransfer.C
    twoResistanceHeatTransfer/twoResistanceHeatTransfer.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/meshTools/lnInclude
    ${LIB_SRC}/finiteVolume/lnInclude
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude
    ${LIB_SRC}/thermophysicalModels/solidThermo/lnInclude
    ${LIB_SRC}/thermophysicalModels/saturationModels/lnInclude
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    ${LIB_SRC}/multiphaseModels/multiphaseProperties/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/solid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/phaseSolid/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    finiteVolume
    sampling
    multicomponentThermophysicalModels
    saturationModels
    compressibleTwoPhases
    multiphaseProperties
    phaseCompressibleMomentumTransportModels
    phaseFluidThermophysicalTransportModels
    phaseFluidMulticomponentThermophysicalTransportModels
    phaseSolidThermophysicalTransportModels
    reactionModels
    fluidThermophysicalTransportModel
    solidThermo
    solidThermophysicalTransportModels
    twoPhaseMixture
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


