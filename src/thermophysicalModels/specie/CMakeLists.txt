lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)
add_library(specie SHARED
    ${sources_from_wmake}
)

target_include_directories(specie PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(specie PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_compile_definitions(specie PRIVATE
    LIB_NAME=libspecie.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};specie" PARENT_SCOPE)


