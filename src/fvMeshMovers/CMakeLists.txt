#add_subdirectory(motionSolver)
add_subdirectory(fvMotionSolvers)
add_subdirectory(interpolator)
add_subdirectory(multiValveEngine)
add_subdirectory(inkJet)

set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)