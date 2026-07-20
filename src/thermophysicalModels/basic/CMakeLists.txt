lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)
add_library(fluidThermophysicalModels SHARED
    ${sources_from_wmake}
    basicThermo/fieldDelegate.h
)
target_include_directories(fluidThermophysicalModels PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(fluidThermophysicalModels PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

)


target_link_libraries(fluidThermophysicalModels PUBLIC
    physicalProperties
    specie
    thermophysicalProperties
    finiteVolume
    meshTools
)

target_compile_definitions(fluidThermophysicalModels PRIVATE
    LIB_NAME=libfluidThermophysicalModels.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};fluidThermophysicalModels" PARENT_SCOPE)


