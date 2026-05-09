lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)
add_library(meshTools SHARED
    ${sources_from_wmake}
)

target_include_directories(meshTools PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(meshTools PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
    triSurface
    fileFormats
    Pstream_best
)


target_compile_definitions(meshTools PRIVATE
    LIB_NAME=libmeshTools.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};meshTools" PARENT_SCOPE)


