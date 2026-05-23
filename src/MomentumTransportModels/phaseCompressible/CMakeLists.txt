lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name phaseCompressibleMomentumTransportModels)

add_library(${target_name} SHARED
    phaseCompressibleMomentumTransportModel.C
    phaseCompressibleMomentumTransportModels.C

    derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonmWallFunction/epsilonmWallFunctionFvPatchScalarField.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
#    ../momentumTransportModels/lnInclude
#    ../compressible/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)

target_link_libraries(${target_name} PUBLIC
    compressibleMomentumTransportModels
#    physicalProperties
#    finiteVolume
#    meshTools
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


