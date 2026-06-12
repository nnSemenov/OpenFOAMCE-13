lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name fluidMulticomponentThermophysicalTransportModels)

add_library(${target_name} SHARED
    fluidMulticomponentThermophysicalTransportModels.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)


target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC

    fluidThermoThermophysicalTransportModels
    fluidThermophysicalModels
    multicomponentThermophysicalModels
    momentumTransportModels
#    specie
#    finiteVolume
#    meshTools
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


