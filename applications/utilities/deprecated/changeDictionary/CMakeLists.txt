lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name changeDictionary)

#wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)

add_executable(${target_name}
#    ${sources_from_wmake}
    changeDictionary.C
)

target_include_directories(${target_name} PRIVATE
    lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OSspecific
    OpenFOAM


    finiteVolume
    polyTopoChange
    meshTools
    Pstream_best
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=libNULL.so
)

set(FOAM_regular_executables "${FOAM_regular_executables};${target_name}" PARENT_SCOPE)
