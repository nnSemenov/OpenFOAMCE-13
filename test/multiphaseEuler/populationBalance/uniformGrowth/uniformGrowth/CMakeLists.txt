cmake_minimum_required(VERSION 3.28)
project(test_find_Mikeno LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Mikeno CONFIG REQUIRED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY $ENV{FOAM_USER_LIBBIN})


add_library(uniformGrowth SHARED
    uniformGrowthGroupFractionFvScalarFieldSource.C
    uniformGrowth.C
)

target_include_directories(uniformGrowth PRIVATE
    $ENV{LIB_SRC}/finiteVolume/lnInclude
    #    ${LIB_SRC}/meshTools/lnInclude
    #    ${LIB_SRC}/sampling/lnInclude
    $ENV{LIB_SRC}/physicalProperties/lnInclude
    $ENV{LIB_SRC}/thermophysicalModels/solidThermo/lnInclude
#    $ENV{LIB_SRC}/thermophysicalModels/specie/lnInclude
    $ENV{LIB_SRC}/thermophysicalModels/basic/lnInclude
    $ENV{LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    $ENV{LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    $ENV{LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    $ENV{LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    #    ${LIB_SRC}/fvModels/general/lnInclude
    #    ${FOAM_MODULES}/multiphaseEuler/phaseSystem/lnInclude
    #    ${FOAM_MODULES}/multiphaseEuler/populationBalance/lnInclude
)

target_link_libraries(uniformGrowth PRIVATE
    Mikeno::OpenFOAM_Defines
    Mikeno::OpenFOAM
    Mikeno::OSspecific

    Mikeno::finiteVolume
    Mikeno::sampling
    Mikeno::meshTools
    Mikeno::momentumTransportModels
    Mikeno::compressibleMomentumTransportModels
    Mikeno::fluidThermoThermophysicalTransportModels
    Mikeno::fvModels
    Mikeno::phaseSystem
    Mikeno::populationBalance
)
