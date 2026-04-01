lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name LagrangianThermo)

add_library(${target_name} SHARED
    basicLagrangianThermo/basicLagrangianThermo.C
    basicLagrangianThermo/BasicLagrangianThermoName.C

    pureLagrangianThermo/pureLagrangianThermo.C

    solidLagrangianThermo/solidLagrangianThermo.C
    solidLagrangianThermo/solidLagrangianThermos.C

    fluidLagrangianThermo/fluidLagrangianThermo.C
    fluidLagrangianThermo/fluidLagrangianThermos.C

    fields/derivedLagrangianFieldSources/density/densityLagrangianScalarFieldSource.C
    fields/derivedLagrangianFieldSources/energy/energyLagrangianScalarFieldSource.C
    fields/derivedLagrangianFieldSources/specificHeatCapacity/specificHeatCapacityLagrangianScalarFieldSource.C
    fields/derivedLagrangianFieldSources/thermalConductivity/thermalConductivityLagrangianScalarFieldSource.C
    fields/derivedLagrangianFieldSources/pressure/pressureLagrangianScalarFieldSource.C
    fields/derivedLagrangianFieldSources/dynamicViscosity/dynamicViscosityLagrangianScalarFieldSource.C
    fields/derivedLagrangianFieldSources/compressibility/compressibilityLagrangianScalarFieldSource.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    physicalProperties
    specie
    fluidThermophysicalModels
    solidThermo
    finiteVolume
    meshTools
    Lagrangian
    Pstream_best
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


