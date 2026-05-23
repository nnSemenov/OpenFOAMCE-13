lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name momentumTransportModels)

add_library(${target_name} SHARED
    momentumTransportModel.C
    RAS/v2f/v2fBase.C
)

set(LESdelta LES/LESdeltas)
target_sources(${target_name} PRIVATE
    ${LESdelta}/LESdelta/LESdelta.C
    ${LESdelta}/cubeRootVolDelta/cubeRootVolDelta.C
    ${LESdelta}/PrandtlDelta/PrandtlDelta.C
    ${LESdelta}/vanDriestDelta/vanDriestDelta.C
    ${LESdelta}/smoothDelta/smoothDelta.C
    ${LESdelta}/maxDeltaxyz/maxDeltaxyz.C
    ${LESdelta}/IDDESDelta/IDDESDelta.C
)

set(LESfilters LES/LESfilters)
target_sources(${target_name} PRIVATE
    ${LESfilters}/LESfilter/LESfilter.C
    ${LESfilters}/simpleFilter/simpleFilter.C
    ${LESfilters}/laplaceFilter/laplaceFilter.C
    ${LESfilters}/anisotropicFilter/anisotropicFilter.C
)

target_sources(${target_name} PRIVATE
    # Turbulence BCs
    derivedFvPatchFields/fixedShearStress/fixedShearStressFvPatchVectorField.C
    derivedFvPatchFields/porousBafflePressure/porousBafflePressureFvPatchField.C
)

# Wall function BCs
set(wallFunctions derivedFvPatchFields/wallFunctions)
target_sources(${target_name} PRIVATE
    ${wallFunctions}/wallCellWallFunction/wallCellWallFunctionFvPatchScalarField.C
)

set(nutWallFunctions ${wallFunctions}/nutWallFunctions)
target_sources(${target_name} PRIVATE
    ${nutWallFunctions}/nutWallFunction/nutWallFunctionFvPatchScalarField.C

    ${nutWallFunctions}/nutkWallFunction/nutkWallFunctionFvPatchScalarField.C
    ${nutWallFunctions}/nutkRoughWallFunction/nutkRoughWallFunctionFvPatchScalarField.C

    ${nutWallFunctions}/nutUWallFunction/nutUWallFunctionFvPatchScalarField.C
    ${nutWallFunctions}/nutUSpaldingWallFunction/nutUSpaldingWallFunctionFvPatchScalarField.C
    ${nutWallFunctions}/nutLowReWallFunction/nutLowReWallFunctionFvPatchScalarField.C
    ${nutWallFunctions}/nutURoughWallFunction/nutURoughWallFunctionFvPatchScalarField.C
)

set(epsilonWallFunctions ${wallFunctions}/epsilonWallFunctions)
target_sources(${target_name} PRIVATE
    ${epsilonWallFunctions}/epsilonWallFunction/epsilonWallFunctionFvPatchScalarField.C
)

set(omegaWallFunctions ${wallFunctions}/omegaWallFunctions)
target_sources(${target_name} PRIVATE
    ${omegaWallFunctions}/omegaWallFunction/omegaWallFunctionFvPatchScalarField.C
)

set(kqRWallFunctions ${wallFunctions}/kqRWallFunctions)
target_sources(${target_name} PRIVATE
    ${kqRWallFunctions}/kqRWallFunction/kqRWallFunctionFvPatchFields.C
    ${kqRWallFunctions}/kLowReWallFunction/kLowReWallFunctionFvPatchScalarField.C
)

set(v2WallFunctions ${wallFunctions}/v2WallFunctions)
target_sources(${target_name} PRIVATE
    ${v2WallFunctions}/v2WallFunction/v2WallFunctionFvPatchScalarField.C
)

set(fWallFunctions ${wallFunctions}/fWallFunctions)
target_sources(${target_name} PRIVATE
    ${fWallFunctions}/fWallFunction/fWallFunctionFvPatchScalarField.C
)

target_sources(${target_name} PRIVATE
    # Inlet turbulence BCs
    RAS/derivedFvPatchFields/turbulentEpsilon/turbulentEpsilonFvPatchScalarField.C
    RAS/derivedFvPatchFields/turbulentOmega/turbulentOmegaFvPatchScalarField.C


    # Turbulence model values
    RAS/derivedFvFieldSources/turbulentEpsilon/turbulentEpsilonFvScalarFieldSource.C
    RAS/derivedFvFieldSources/turbulentOmega/turbulentOmegaFvScalarFieldSource.C
)

set(generalisedNewtonianViscosityModels laminar/generalisedNewtonian/generalisedNewtonianViscosityModels)
target_sources(${target_name} PRIVATE
    ${generalisedNewtonianViscosityModels}/generalisedNewtonianViscosityModel/generalisedNewtonianViscosityModel.C
    ${generalisedNewtonianViscosityModels}/generalisedNewtonianViscosityModel/generalisedNewtonianViscosityModelNew.C

    ${generalisedNewtonianViscosityModels}/Newtonian/NewtonianViscosityModel.C
)

set(strainRateViscosityModels ${generalisedNewtonianViscosityModels}/strainRateViscosityModels)
target_sources(${target_name} PRIVATE
    ${strainRateViscosityModels}/strainRateViscosityModel/strainRateViscosityModel.C
    ${strainRateViscosityModels}/CrossPowerLaw/CrossPowerLaw.C
    ${strainRateViscosityModels}/BirdCarreau/BirdCarreau.C
    ${strainRateViscosityModels}/Casson/Casson.C
    ${strainRateViscosityModels}/HerschelBulkley/HerschelBulkley.C
    ${strainRateViscosityModels}/powerLaw/powerLaw.C
    ${strainRateViscosityModels}/strainRateFunction/strainRateFunction.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)

target_link_libraries(${target_name} PUBLIC
    finiteVolume
    meshTools
    Pstream_best
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


