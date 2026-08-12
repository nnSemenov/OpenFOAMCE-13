lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name porousMediaFluidSolver)

add_library(${target_name} SHARED
    porousMediaFluid.C
    porousPhaseInfo.C
    thermophysicalPredictor.C
    effectiveHeatConductivity.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_compile_features(${target_name} PUBLIC
    cxx_std_23
)

target_include_directories(${target_name} PRIVATE
    #    ${LIB_SRC}/finiteVolume/lnInclude
    #    ${LIB_SRC}/meshTools/lnInclude
    #    ${LIB_SRC}/sampling/lnInclude
    #    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluidThermo/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/solid/lnInclude
    #    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    Boost::boost
)

target_link_libraries(${target_name} PRIVATE
    fluidSolver

    specie
    solidThermo
    solidThermophysicalTransportModels
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


