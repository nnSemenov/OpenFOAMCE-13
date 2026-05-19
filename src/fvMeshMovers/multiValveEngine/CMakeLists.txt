lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name fvMeshMoversMultiValveEngine)

add_library(${target_name} SHARED
    multiValveEngine.C
    movingObject.C
    piston.C
    valve.C
    valveList.C
    crankConnectingRodMotion/crankConnectingRodMotion.C
    pistonPointEdgeData/pistonPointEdgeData.C
    multiValveEngineState/multiValveEngineState.C
    zoneGenerators/cylinderHeadPoints/cylinderHeadPoints.C
    zoneGenerators/pistonBowlPoints/pistonBowlPoints.C
    solidBodyMotionFunctions/multiValveEnginePistonMotion/multiValveEnginePistonMotion.C
    solidBodyMotionFunctions/multiValveEngineValveMotion/multiValveEngineValveMotion.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC
    finiteVolume
    meshTools
    pointMeshMovers
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


