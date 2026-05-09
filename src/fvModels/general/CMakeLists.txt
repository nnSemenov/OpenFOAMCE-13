lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name fvModels)

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)

add_library(${target_name} SHARED
    ${sources_from_wmake}
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)


target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    finiteVolume 
    sampling 
    meshTools 
    momentumTransportModels 
    compressibleMomentumTransportModels
    fluidThermophysicalTransportModel
    fluidThermoThermophysicalTransportModels 
    physicalProperties 
    solidThermo 
    multicomponentThermophysicalModels 
    fluidMulticomponentThermophysicalTransportModels
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


