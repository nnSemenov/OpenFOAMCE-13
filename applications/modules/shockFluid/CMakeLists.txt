lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name shockFluidSolver)

add_library(${target_name} SHARED
    setRDeltaT.C
    moveMesh.C
    fluxPredictor.C
    correctDensity.C
    momentumPredictor.C
    thermophysicalPredictor.C
    pressureCorrector.C
    shockFluid.C

    derivedFvPatchFields/mixedFixedValueSlip/mixedFixedValueSlipFvPatchFields.C
    derivedFvPatchFields/U/maxwellSlipUFvPatchVectorField.C
    derivedFvPatchFields/T/smoluchowskiJumpTFvPatchScalarField.C
    derivedFvPatchFields/rho/fixedRhoFvPatchScalarField.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluidThermo/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    basicFluidSolver

    fluidThermophysicalModels
    fluidThermoThermophysicalTransportModels
    momentumTransportModels
    compressibleMomentumTransportModels
    coupledThermophysicalTransportModels
    finiteVolume
    meshTools
    fvModels
    fvConstraints
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


