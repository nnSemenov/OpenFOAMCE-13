lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(genericFields SHARED
    genericFieldBase/genericFieldBase.C
)

target_include_directories(genericFields PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(genericFields PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    genericPatches
)

target_compile_definitions(genericFields PRIVATE
    LIB_NAME=libgenericFields.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};genericFields" PARENT_SCOPE)


