lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name snappyHexMesh_exe)

add_executable(${target_name}
    snappyHexMesh.C
)

set_target_properties(${target_name} PROPERTIES
    OUTPUT_NAME snappyHexMesh
)

target_include_directories(${target_name} PRIVATE
    lnInclude
    ${LIB_SRC}/triSurface/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OSspecific
    OpenFOAM

    finiteVolume 
    genericFvFields 
    decompositionMethods
    meshTools 
    meshCheck 
    surfMesh 
    fileFormats 
    polyTopoChange 
    snappyHexMesh
    triSurface
    Pstream_best
)

set(extra_deps scotchDecomp ptscotchDecomp)
foreach (dep ${extra_deps})
    if(TARGET ${dep})
        target_link_libraries(${target_name} PRIVATE
            ${dep}
        )
    endif ()
endforeach ()

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=libNULL.so
)

set(FOAM_regular_executables "${FOAM_regular_executables};${target_name}" PARENT_SCOPE)