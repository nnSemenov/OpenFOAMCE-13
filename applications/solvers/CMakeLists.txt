set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${LTO_option})

add_subdirectory(foamRun)
add_subdirectory(foamMultiRun)
add_subdirectory(chemFoam)
add_subdirectory(potentialFoam)
add_subdirectory(boundaryFoam)

set(FOAM_regular_executables ${FOAM_regular_executables} PARENT_SCOPE)
set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)
set(FOAM_special_libraries ${FOAM_special_libraries} PARENT_SCOPE)