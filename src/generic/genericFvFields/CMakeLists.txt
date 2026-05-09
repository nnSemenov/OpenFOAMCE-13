lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(genericFvFields SHARED
    genericFvPatchField/genericFvPatchFields.C
    genericPointPatchField/genericPointPatchFields.C
    genericFvFieldSource/genericFvFieldSources.C
)

target_include_directories(genericFvFields PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(genericFvFields PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    genericPatches
    genericFields
    genericFvPatches
    finiteVolume
)

target_compile_definitions(genericFvFields PRIVATE
    LIB_NAME=libgenericFvFields.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};genericFvFields" PARENT_SCOPE)


