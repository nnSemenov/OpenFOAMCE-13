set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${LTO_option})

add_subdirectory(mesh)
add_subdirectory(miscellaneous)
add_subdirectory(parallelProcessing)
add_subdirectory(postProcessing)
add_subdirectory(preProcessing)
add_subdirectory(surface)
add_subdirectory(thermophysical)

set(FOAM_regular_executables ${FOAM_regular_executables} PARENT_SCOPE)
set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)
set(FOAM_special_libraries "${FOAM_special_libraries}" PARENT_SCOPE)
