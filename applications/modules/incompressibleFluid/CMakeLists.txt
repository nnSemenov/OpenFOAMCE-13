lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name incompressibleFluidSolver)

add_library(${target_name} SHARED
    setRDeltaT.C
    moveMesh.C
    momentumPredictor.C
    correctPressure.C
    incompressibleFluid.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE

    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude 
    ${LIB_SRC}/MomentumTransportModels/incompressible/lnInclude 
    ${LIB_SRC}/finiteVolume/lnInclude 
    ${LIB_SRC}/meshTools/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
    
    basicFluidSolver
    physicalProperties
    incompressibleMomentumTransportModels
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


