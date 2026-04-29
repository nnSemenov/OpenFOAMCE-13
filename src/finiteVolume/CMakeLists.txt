lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)
add_library(finiteVolume SHARED
    ${sources_from_wmake}
)


target_include_directories(finiteVolume PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(finiteVolume PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
    triSurface
    meshTools
    Pstream_best
)


target_compile_definitions(finiteVolume PRIVATE
    LIB_NAME=libfiniteVolume.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};finiteVolume" PARENT_SCOPE)


