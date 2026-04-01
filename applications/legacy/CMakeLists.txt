set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${LTO_option})

add_subdirectory(basic/laplacianFoam)
add_subdirectory(basic/financialFoam)

add_subdirectory(combustion/PDRFoam)

add_subdirectory(compressible/rhoPorousSimpleFoam)

add_subdirectory(incompressible/adjointShapeOptimisationFoam)
add_subdirectory(incompressible/icoFoam)
add_subdirectory(incompressible/porousSimpleFoam)
add_subdirectory(incompressible/shallowWaterFoam)

add_subdirectory(electromagnetics/electrostaticFoam)
add_subdirectory(electromagnetics/magneticFoam)
add_subdirectory(electromagnetics/mhdFoam)

add_subdirectory(lagrangian/dsmcFoam)
add_subdirectory(lagrangian/mdEquilibrationFoam)
add_subdirectory(lagrangian/mdFoam)

set(FOAM_regular_executables ${FOAM_regular_executables} PARENT_SCOPE)
set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)
set(FOAM_special_libraries ${FOAM_special_libraries} PARENT_SCOPE)