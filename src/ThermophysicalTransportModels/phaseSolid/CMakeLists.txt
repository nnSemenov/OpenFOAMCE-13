lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name phaseSolidThermophysicalTransportModels)

add_library(${target_name} SHARED
    phaseSolidThermophysicalTransportModel/phaseSolidThermophysicalTransportModel.C
    phaseSolidThermophysicalTransportModels.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)


target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/solid/lnInclude
    ${LIB_SRC}
    ${LIB_SRC}
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)

target_link_libraries(${target_name} PUBLIC
    thermophysicalTransportModel
    solidThermo
#    specie
#    finiteVolume
#    meshTools
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


