

lnInclude(OpenFOAM RESULT_DIR OpenFOAM_lnInclude_dir)
message(STATUS "OpenFOAM_lnInclude_dir = ${OpenFOAM_lnInclude_dir}")

add_subdirectory(wmakeParse)
add_subdirectory(OSspecific/POSIX)
add_subdirectory(Pstream)
add_subdirectory(OpenFOAM)

add_subdirectory(fileFormats)
add_subdirectory(surfMesh)
add_subdirectory(triSurface)
add_subdirectory(meshTools)

add_subdirectory(finiteVolume)

add_subdirectory(ODE)
add_subdirectory(physicalProperties)
add_subdirectory(thermophysicalModels)

add_subdirectory(tracking)
add_subdirectory(lagrangian)
add_subdirectory(Lagrangian)

add_subdirectory(generic)

add_subdirectory(sampling)
add_subdirectory(meshCheck)
#add_subdirectory(topoSetSources)
#add_subdirectory(fvTopoSetSources)

add_subdirectory(motionSolvers)

add_subdirectory(mesh)
add_subdirectory(polyTopoChange)

add_subdirectory(conversion)

add_subdirectory(parallel)

add_subdirectory(fvMeshStitchers)
add_subdirectory(fvMeshMovers)
add_subdirectory(fvMeshTopoChangers)
add_subdirectory(fvMeshDistributors)

add_subdirectory(randomProcesses)

add_subdirectory(twoPhaseModels)
add_subdirectory(multiphaseModels)
add_subdirectory(MomentumTransportModels)
add_subdirectory(ThermophysicalTransportModels)
add_subdirectory(radiationModels)
add_subdirectory(reactionModels)


add_subdirectory(renumber)
add_subdirectory(fvAgglomerationMethods)
add_subdirectory(fvMotionSolver)

add_subdirectory(fvConstraints)

add_subdirectory(fvModels)
add_subdirectory(functionObjects)

add_subdirectory(rigidBodyMotion)
add_subdirectory(specieTransfer)
add_subdirectory(atmosphericModels)
add_subdirectory(waves)
###############################################

set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)
set(FOAM_special_libraries ${FOAM_special_libraries} PARENT_SCOPE)