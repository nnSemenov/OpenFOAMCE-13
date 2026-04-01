lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name fvModels)

add_library(${target_name} SHARED
    codedFvModel/codedFvModel.C
    semiImplicitSource/semiImplicitSource.C

    actuationDisk/actuationDisk.C
    effectivenessHeatExchanger/effectivenessHeatExchanger.C
    porosityForce/porosityForce.C
    radialActuationDisk/radialActuationDisk.C
    solidificationMelting/solidificationMelting.C
    sixDoFAcceleration/sixDoFAcceleration.C
    buoyancyForce/buoyancyForce.C
    buoyancyEnergy/buoyancyEnergy.C
    viscousHeating/viscousHeating.C
    phaseLimitStabilisation/phaseLimitStabilisation.C
    acceleration/acceleration.C
    volumeBlockage/volumeBlockage.C
    solidThermalEquilibrium/solidThermalEquilibrium.C
    volumeSource/volumeSource.C
    massSource/massSourceBase.C
    massSource/massSource.C
    heatSource/heatSource.C
    zeroDimensionalMassSource/zeroDimensionalMassSourceBase.C
    zeroDimensionalMassSource/zeroDimensionalMassSource.C
    massTransfer/massTransfer.C
    massTransfer/coefficientMassTransfer.C
    phaseChange/phaseChange.C
    phaseChange/coefficientPhaseChange.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)


target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    finiteVolume 
    sampling 
    meshTools 
    momentumTransportModels 
    compressibleMomentumTransportModels
    fluidThermophysicalTransportModel
    fluidThermoThermophysicalTransportModels 
    physicalProperties 
    solidThermo 
    multicomponentThermophysicalModels 
    fluidMulticomponentThermophysicalTransportModels
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


