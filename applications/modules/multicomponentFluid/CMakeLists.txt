lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name multicomponentFluidSolver)

add_library(${target_name} SHARED
    multicomponentFluid.C
    setRDeltaT.C
    thermophysicalPredictor.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${LIB_SRC}/physicalProperties/lnInclude
    #${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC

    isothermalFluidSolver

    fluidThermophysicalModels
    fluidThermophysicalTransportModel
    specie
    chemistryModel
    ODE
    reactionModels
    multicomponentThermophysicalModels
    momentumTransportModels
    compressibleMomentumTransportModels
    coupledThermophysicalTransportModels
    fluidMulticomponentThermophysicalTransportModels
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


