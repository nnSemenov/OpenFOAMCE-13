lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(genericPatches SHARED
    genericPolyPatch/genericPolyPatch.C
)

target_include_directories(genericPatches PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(genericPatches PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_compile_definitions(genericPatches PRIVATE
    LIB_NAME=libgenericPatches.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};genericPatches" PARENT_SCOPE)


