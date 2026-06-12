lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name phaseIncompressibleMomentumTransportModels)

add_library(${target_name} SHARED
    phaseIncompressibleMomentumTransportModel.C
    phaseIncompressibleMomentumTransportModels.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
#    ../momentumTransportModels/lnInclude
#    ../incompressible/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

target_link_libraries(${target_name} PUBLIC
    incompressibleMomentumTransportModels
#    physicalProperties
#    finiteVolume
#    meshTools
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


