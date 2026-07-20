lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)
add_library(thermophysicalProperties SHARED
    ${sources_from_wmake}
)
target_include_directories(thermophysicalProperties PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(thermophysicalProperties PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    specie
)

target_compile_definitions(thermophysicalProperties PRIVATE
    LIB_NAME=libthermophysicalProperties.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};thermophysicalProperties" PARENT_SCOPE)


