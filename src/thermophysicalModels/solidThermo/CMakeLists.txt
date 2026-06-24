lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)

add_library(solidThermo SHARED
    ${sources_from_wmake}
)
target_include_directories(solidThermo PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(solidThermo PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
)


target_link_libraries(solidThermo PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    meshTools 
    finiteVolume 
    specie 
    fluidThermophysicalModels 
    multicomponentThermophysicalModels
)

target_compile_definitions(solidThermo PRIVATE
    LIB_NAME=libsolidThermo.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};solidThermo" PARENT_SCOPE)


