lnInclude(${CMAKE_CURRENT_SOURCE_DIR} HEADER_LIST target_headers)

set(target_name phaseSystem)

wmake_parse_files(${CMAKE_CURRENT_SOURCE_DIR}/Make/files sources_from_wmake)
add_library(${target_name} SHARED
    ${target_headers}

    ${sources_from_wmake}
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    #    ${LIB_SRC}/meshTools/lnInclude
    #    ${LIB_SRC}/finiteVolume/lnInclude
    #    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude
    ${LIB_SRC}/thermophysicalModels/solidThermo/lnInclude
    ${LIB_SRC}/thermophysicalModels/saturationModels/lnInclude
    #    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    #    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    #    ${LIB_SRC}/multiphaseModels/multiphaseProperties/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/solid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/phaseSolid/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_link_libraries(${target_name} PUBLIC
    finiteVolume
    meshTools
    sampling
    physicalProperties
    multicomponentThermophysicalModels
    saturationModels
    compressibleTwoPhases
    multiphaseProperties
    phaseCompressibleMomentumTransportModels
    phaseFluidThermophysicalTransportModels
    phaseFluidMulticomponentThermophysicalTransportModels
    phaseSolidThermophysicalTransportModels
    reactionModels
    fluidThermophysicalTransportModel
    solidThermo
    solidThermophysicalTransportModels
    twoPhaseMixture
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


