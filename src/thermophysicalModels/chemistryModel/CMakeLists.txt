lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)

add_library(chemistryModel SHARED
    ${sources_from_wmake}
)

target_include_directories(chemistryModel PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(chemistryModel PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/functions/Polynomial
)

target_link_libraries(chemistryModel PUBLIC
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    fluidThermophysicalModels 
    multicomponentThermophysicalModels 
    specie 
    ODE 
    finiteVolume 
    meshTools
)

target_compile_definitions(chemistryModel PRIVATE
    LIB_NAME=libchemistryModel.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};chemistryModel" PARENT_SCOPE)


