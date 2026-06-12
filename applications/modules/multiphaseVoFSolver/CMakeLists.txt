lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name multiphaseVoFSolver)

add_library(${target_name} SHARED
    multiphaseVoFMixture/VoFphase/VoFphase.C
    multiphaseVoFMixture/multiphaseVoFMixture.C
    setInterfaceRDeltaT.C
    #alphaPredictor.C
    multiphaseVoFSolver.C
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
)

target_link_libraries(${target_name} PUBLIC
    VoFSolver
    #    basicFluidSolver
    #    physicalProperties
    multiphaseProperties
    #    interfaceCompression
    #    finiteVolume
    #    meshTools
    #    fvModels
    #    fvConstraints
    #    sampling
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


