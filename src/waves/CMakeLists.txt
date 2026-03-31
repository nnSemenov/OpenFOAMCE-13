lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name waves)

add_library(${target_name} SHARED
    waveModels/waveModel/waveModel.C
    waveModels/waveModel/waveModelNew.C
    waveModels/Airy/AiryCoeffs.C
    waveModels/Airy/Airy.C
    waveModels/Stokes2/Stokes2.C
    waveModels/Stokes5/Stokes5.C
    waveModels/solitary/solitary.C
    waveModels/irregular/waveSpectra/waveSpectrum/waveSpectrum.C
    waveModels/irregular/waveSpectra/PiersonMoskowitz/PiersonMoskowitz.C
    waveModels/irregular/waveSpectra/JONSWAP/JONSWAP.C
    waveModels/irregular/waveSpectra/GodaJONSWAP/GodaJONSWAP.C
    waveModels/irregular/irregular.C

    waveSuperpositions/waveSuperposition/waveSuperposition.C
    waveSuperpositions/waveSuperposition/waveSuperpositionNew.C

    dimensionedFieldFunctions/waveAlpha/waveAlpha_DimensionedFieldFunction.C
    dimensionedFieldFunctions/waveVelocity/waveVelocity_DimensionedFieldFunction.C

    derivedFvPatchFields/waveAlpha/waveAlphaFvPatchScalarField.C
    derivedFvPatchFields/waveInletOutlet/waveInletOutletFvPatchFields.C
    derivedFvPatchFields/waveVelocity/waveVelocityFvPatchVectorField.C

    fvModels/forcing/forcing.C
    fvModels/isotropicDamping/isotropicDamping.C
    fvModels/verticalDamping/verticalDamping.C
    fvModels/waveForcing/waveForcing.C
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

    finiteVolume 
    meshTools 
    polyTopoChange 
    atmosphericModels 
    sampling
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


