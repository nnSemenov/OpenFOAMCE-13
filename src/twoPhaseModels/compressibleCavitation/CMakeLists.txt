lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name compressibleCavitationModels)

add_library(${target_name} SHARED
    cavitationModel/cavitationModel.C
    cavitationModel/cavitationModelNew.C
    Kunz/Kunz.C
    Merkle/Merkle.C
    SchnerrSauer/SchnerrSauer.C
    Saito/Saito.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    compressibleTwoPhases 
    physicalProperties 
    fluidThermophysicalModels
    finiteVolume 
    meshTools
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


