lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(fluidThermophysicalModels SHARED
    basicThermo/basicThermo.C
    basicThermo/BasicThermoName.C

    pureThermo/pureThermo.C

    rhoThermo/rhoThermo.C

    fluidThermo/fluidThermo.C
    fluidThermo/hydrostaticInitialisation.C

    psiThermo/psiThermo.C
    psiThermo/psiThermos.C

    rhoFluidThermo/rhoFluidThermo.C
    rhoFluidThermo/rhoFluidThermos.C

    liquidThermo/liquidThermo.C
    liquidThermo/liquidThermos.C

    derivedFvPatchFields/energyJump/energyJumpFvPatchScalarField.C
    derivedFvPatchFields/fixedEnergy/fixedEnergyFvPatchScalarField.C
    derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.C
    derivedFvPatchFields/gradientEnergy/gradientEnergyCalculatedTemperatureFvPatchScalarField.C
    derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.C
    derivedFvPatchFields/mixedEnergy/mixedEnergyCalculatedTemperatureFvPatchScalarField.C

    derivedFvFieldSources/energy/energyFvScalarFieldSource.C
    derivedFvFieldSources/energy/energyCalculatedTemperatureFvScalarFieldSource.C
    derivedFvFieldSources/uniformFixedEnergyTemperature/uniformFixedEnergyTemperatureFvScalarFieldSource.C
    derivedFvFieldSources/uniformInletOutletEnergyTemperature/uniformInletOutletEnergyTemperatureFvScalarFieldSource.C

    basicThermo/fieldDelegate.h
)
target_include_directories(fluidThermophysicalModels PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(fluidThermophysicalModels PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)


target_link_libraries(fluidThermophysicalModels PUBLIC
    physicalProperties
    specie
    thermophysicalProperties
    finiteVolume
    meshTools
)

target_compile_definitions(fluidThermophysicalModels PRIVATE
    LIB_NAME=libfluidThermophysicalModels.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};fluidThermophysicalModels" PARENT_SCOPE)


