lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name atmosphericModels)

add_library(${target_name} SHARED
    atmosphericBoundaryLayer/atmosphericBoundaryLayer.C

    derivedFvPatchFields/atmosphericBoundaryLayerVelocity/atmosphericBoundaryLayerVelocityFvPatchVectorField.C
    derivedFvPatchFields/atmosphericBoundaryLayerTurbulentKineticEnergy/atmosphericBoundaryLayerTurbulentKineticEnergyFvPatchScalarField.C
    derivedFvPatchFields/atmosphericBoundaryLayerTurbulentEpsilon/atmosphericBoundaryLayerTurbulentEpsilonFvPatchScalarField.C
    derivedFvPatchFields/nutAtmosphericBoundaryLayerWallFunction/nutAtmosphericBoundaryLayerWallFunctionFvPatchScalarField.C

    DimensionedFieldFunctions/makeDimensionedFieldFunctions.C

    atmosphericTurbulentTransportModels.C
    porosityModels/powerLawLopesdaCosta/powerLawLopesdaCosta.C
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

    momentumTransportModels 
    incompressibleMomentumTransportModels 
    physicalProperties 
    finiteVolume 
    meshTools 
    sampling 
    triSurface 
    fvModels 
    fvConstraints
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


