lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name forces)

add_library(${target_name} SHARED
    forcesBase/forcesBase.C
    forces/forces.C
    movingForces/movingForces.C
    forceCoeffs/forceCoeffs.C

    sectionalForcesBase/sectionalForcesBase.C
    sectionalForceGraph/sectionalForceGraph.C
    sectionalForceProbes/sectionalForceProbes.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE

)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    fluidThermophysicalModels 
    physicalProperties 
    momentumTransportModels 
    incompressibleMomentumTransportModels 
    compressibleMomentumTransportModels 
    phaseIncompressibleMomentumTransportModels 
    phaseCompressibleMomentumTransportModels 
    specie 
    fileFormats
    sampling
    finiteVolume 
    meshTools
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


