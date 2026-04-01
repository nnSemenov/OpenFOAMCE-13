add_subdirectory(Lagrangian)
add_subdirectory(LagrangianThermo)
add_subdirectory(LagrangianFunctionObjects)
add_subdirectory(cloud)
add_subdirectory(cloudFunctionObjects)

set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)

