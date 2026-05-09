lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(genericLagrangianFields SHARED
    genericLagrangianPatchField/genericLagrangianPatchFields.C
    genericLagrangianFieldSource/genericLagrangianFieldSources.C
)

target_include_directories(genericLagrangianFields PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(genericLagrangianFields PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    genericPatches
    genericFields
    genericLagrangianPatches
    Lagrangian
    finiteVolume
)

target_compile_definitions(genericLagrangianFields PRIVATE
    LIB_NAME=libgenericLagrangianFields.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};genericLagrangianFields" PARENT_SCOPE)


