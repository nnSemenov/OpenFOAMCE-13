lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(genericFvPatches SHARED
    genericFvPatch/genericFvPatch.C
    genericPointPatch/genericPointPatch.C
)

target_include_directories(genericFvPatches PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(genericFvPatches PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    genericPatches
    finiteVolume
)

target_compile_definitions(genericFvPatches PRIVATE
    LIB_NAME=libgenericFvPatches.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};genericFvPatches" PARENT_SCOPE)


