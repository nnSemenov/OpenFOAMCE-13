lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(saturationModels SHARED
    saturationPressureModel/saturationPressureModel.C
    saturationPressureModel/saturationPressureModelNew.C
    saturationTemperatureModel/saturationTemperatureModel.C
    saturationTemperatureModel/saturationTemperatureModelNew.C

    Antoine/Antoine.C
    AntoineExtended/AntoineExtended.C
    ArdenBuck/ArdenBuck.C
    constantPressure/constantPressure.C
    constantTemperature/constantTemperature.C
    function1Temperature/function1Temperature.C
    function1Pressure/function1Pressure.C
    polynomialTemperature/polynomialTemperature.C
)

target_include_directories(saturationModels PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(saturationModels PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/thermophysicalProperties/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude 
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude 
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude 
    ${LIB_SRC}/finiteVolume/lnInclude 
    ${LIB_SRC}/meshTools/lnInclude 
    ${LIB_SRC}/sampling/lnInclude
)

target_link_libraries(saturationModels PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    fluidThermophysicalModels 
    multicomponentThermophysicalModels 
    specie
)

target_compile_definitions(saturationModels PRIVATE
    LIB_NAME=libsaturationModels.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};saturationModels" PARENT_SCOPE)


