if(NOT ${MIKENO_UNIT_TEST})
    message(STATUS "Unit tests are skipped.")
    return()
endif ()

enable_testing()

function(FOAM_unit_test)
    cmake_parse_arguments(arg "BUILD_ONLY;WILL_FAIL" "NAME;TEST_NAME;WORKING_DIRECTORY" "SOURCES;LINK_LIBS;INCLUDE_DIRS;TEST_COMMAND_ARGS" ${ARGN})

    list(GET arg_SOURCES 0 first_source_file)
    cmake_path(GET first_source_file PARENT_PATH test_proj_dir)
    if(NOT DEFINED arg_NAME)
        cmake_path(GET test_proj_dir FILENAME arg_NAME)
        set(arg_NAME Test-${arg_NAME})
    endif ()
    add_executable(${arg_NAME} ${arg_SOURCES})
    target_include_directories(${arg_NAME} PRIVATE
        ${arg_INCLUDE_DIRS}
    )
    target_link_libraries(${arg_NAME} PRIVATE
        OpenFOAM_Defines
        OpenFOAM
        OSspecific

        ${arg_LINK_LIBS}
    )

    if(${arg_BUILD_ONLY})
        return()
    endif ()

    if(NOT DEFINED arg_TEST_NAME)
        set(arg_TEST_NAME ${arg_NAME})
    endif ()

    cmake_path(ABSOLUTE_PATH test_proj_dir)
    if(NOT DEFINED arg_WORKING_DIRECTORY)
        set(arg_WORKING_DIRECTORY ${test_proj_dir})
    endif ()

    add_test(NAME ${arg_TEST_NAME}
        COMMAND ${arg_NAME} ${arg_TEST_COMMAND_ARGS}
        WORKING_DIRECTORY ${arg_WORKING_DIRECTORY}
    )
    if(${arg_WILL_FAIL})
        set_property(TEST ${arg_TEST_NAME}
            PROPERTY WILL_FAIL TRUE
        )
    endif ()
endfunction()

FOAM_unit_test(SOURCES alloc/Test-alloc.C)
FOAM_unit_test(SOURCES BinSum/Test-BinSum.C LINK_LIBS Pstream_best)
FOAM_unit_test(SOURCES boundSphere/Test-boundSphere.C
    LINK_LIBS fileFormats Pstream_best
    TEST_COMMAND_ARGS -nTests 1000 1000
)
FOAM_unit_test(SOURCES callback/Test-callback.C)
FOAM_unit_test(SOURCES Circulator/Test-Circulator.C)
FOAM_unit_test(SOURCES codeStream/Test-codeStream.C
    TEST_COMMAND_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/codeStream/codeStreamDict1
)
FOAM_unit_test(SOURCES CompactIOList/Test-CompactIOList.C
    TEST_COMMAND_ARGS -case ${PROJECT_SOURCE_DIR}/tutorials/incompressibleFluid/mixerSRF
)
FOAM_unit_test(SOURCES CompactListList/Test-CompactListList.C
    TEST_COMMAND_ARGS -case ${PROJECT_SOURCE_DIR}/tutorials/incompressibleFluid/mixerSRF
)

FOAM_unit_test(SOURCES cubicEqn/Test-cubicEqn.C)

FOAM_unit_test(SOURCES cyclic/Test-cyclic.C
    LINK_LIBS finiteVolume meshTools
    TEST_COMMAND_ARGS -case ${PROJECT_SOURCE_DIR}/tutorials/incompressibleFluid/mixerSRF
)

FOAM_unit_test(SOURCES decomposedBlockData/Test-decomposedBlockData.C
    BUILD_ONLY
)

FOAM_unit_test(SOURCES delete/Test-delete.C)

FOAM_unit_test(SOURCES Dictionary/Test-Dictionary.C)

FOAM_unit_test(SOURCES dictionary/Test-dictionary.C
    BUILD_ONLY
)

FOAM_unit_test(SOURCES dimensionedType/Test-dimensionedType.C
    BUILD_ONLY
)

FOAM_unit_test(SOURCES distribution/Test-distribution.C
    LINK_LIBS sampling
    BUILD_ONLY
)

FOAM_unit_test(SOURCES DLList/Test-DLList.C)

FOAM_unit_test(SOURCES DynamicField/Test-DynamicField.C)

FOAM_unit_test(SOURCES dynamicIndexedOctree/Test-dynamicIndexedOctree.C
    LINK_LIBS meshTools
)

FOAM_unit_test(SOURCES DynamicList/Test-DynamicList.C)

FOAM_unit_test(SOURCES error/Test-error.C
#    WILL_FAIL
    BUILD_ONLY
)

FOAM_unit_test(SOURCES extendedStencil/Test-ExtendedStencil.C
    LINK_LIBS finiteVolume meshTools Pstream_best
    TEST_COMMAND_ARGS -case ${PROJECT_SOURCE_DIR}/tutorials/mesh/blockMesh/pipe
)

FOAM_unit_test(SOURCES fieldDependency/Test-fieldDependency.C
    LINK_LIBS finiteVolume
    TEST_COMMAND_ARGS -case ${PROJECT_SOURCE_DIR}/tutorials/incompressibleFluid/mixerSRF
)

FOAM_unit_test(SOURCES fileName/Test-fileName.C)

FOAM_unit_test(SOURCES fileNameClean/Test-fileNameClean.C)













FOAM_unit_test(SOURCES realFluidEOS/Test-realFluidEOS.cpp
    LINK_LIBS specie
)
