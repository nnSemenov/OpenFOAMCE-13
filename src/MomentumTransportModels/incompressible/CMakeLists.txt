lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name incompressibleMomentumTransportModels)

add_library(${target_name} SHARED
    incompressibleMomentumTransportModel.C
    incompressibleMomentumTransportModels.C

    RAS/qZeta/qZeta.C
    RAS/kkLOmega/kkLOmega.C
    RAS/LamBremhorstKE/LamBremhorstKE.C
    RAS/LienLeschziner/LienLeschziner.C
    RAS/ShihQuadraticKE/ShihQuadraticKE.C
    RAS/LienCubicKE/LienCubicKE.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)

target_link_libraries(${target_name} PUBLIC
    physicalProperties
    momentumTransportModels
#    finiteVolume
#    meshTools
#    Pstream_best
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


