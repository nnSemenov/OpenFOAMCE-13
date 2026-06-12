

set(processed_file ${CMAKE_CURRENT_BINARY_DIR}/global/global.C)
process_Cver(global/global.Cver ${processed_file})
add_library(OpenFOAM SHARED
    ${processed_file}
)

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files openfoam_sources)
target_sources(OpenFOAM PRIVATE
    ${openfoam_sources}
)
#message(STATUS "Openfoam sources: ${openfoam_sources}")

if (NOT ${WM_PRECISION_OPTION} STREQUAL "SP")
    target_sources(OpenFOAM PRIVATE
        primitives/Vector/floatVector/floatVector.C
        primitives/Tensor/floatTensor/floatTensor.C
    )
endif ()


target_link_libraries(OpenFOAM PRIVATE
    OSspecific
    wmakeParse::wmakeParse
    Pstream_best
    ZLIB::ZLIB
)

target_compile_definitions(OpenFOAM PRIVATE
    LIB_NAME=libOpenFOAM.so
    WM_OPTIONS=${WM_OPTIONS}
)

if (${LINUX})
    target_link_options(OpenFOAM PRIVATE
        "-Wl,--no-allow-shlib-undefined"
    )
endif ()

target_include_directories(OpenFOAM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude
)


set(FOAM_regular_libraries "${FOAM_regular_libraries};OpenFOAM" PARENT_SCOPE)
