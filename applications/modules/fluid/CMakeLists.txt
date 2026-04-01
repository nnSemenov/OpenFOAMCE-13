lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name fluidSolver)

add_library(${target_name} SHARED
    fluid.C
    thermophysicalPredictor.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/finiteVolume/cfdTools
    ${LIB_SRC}/finiteVolume/lnInclude
    ${LIB_SRC}/meshTools/lnInclude
    ${LIB_SRC}/sampling/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluidThermo/lnInclude

)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    basicFluidSolver
    isothermalFluidSolver

    fluidThermophysicalModels 
    thermophysicalTransportModel
    fluidThermophysicalTransportModel
    fluidThermoThermophysicalTransportModels
    compressibleMomentumTransportModels 
    coupledThermophysicalTransportModels
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


