lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name multiphaseEulerMomentumTransportModels)

add_library(${target_name} SHARED
    momentumTransportModels.C

    phasePressureModel/phasePressureModel.C

    kineticTheoryModels/kineticTheoryModel/kineticTheoryModel.C

    kineticTheoryModels/viscosityModel/viscosityModel/kineticTheoryViscosityModel.C
    kineticTheoryModels/viscosityModel/viscosityModel/kineticTheoryViscosityModelNew.C
    kineticTheoryModels/viscosityModel/Gidaspow/GidaspowViscosity.C
    kineticTheoryModels/viscosityModel/Syamlal/SyamlalViscosity.C
    kineticTheoryModels/viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C
    kineticTheoryModels/viscosityModel/none/noneViscosity.C

    kineticTheoryModels/conductivityModel/conductivityModel/conductivityModel.C
    kineticTheoryModels/conductivityModel/conductivityModel/conductivityModelNew.C
    kineticTheoryModels/conductivityModel/Gidaspow/GidaspowConductivity.C
    kineticTheoryModels/conductivityModel/Syamlal/SyamlalConductivity.C
    kineticTheoryModels/conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C

    kineticTheoryModels/radialModel/radialModel/radialModel.C
    kineticTheoryModels/radialModel/radialModel/radialModelNew.C
    kineticTheoryModels/radialModel/CarnahanStarling/CarnahanStarlingRadial.C
    kineticTheoryModels/radialModel/LunSavage/LunSavageRadial.C
    kineticTheoryModels/radialModel/SinclairJackson/SinclairJacksonRadial.C

    kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModel.C
    kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModelNew.C
    kineticTheoryModels/granularPressureModel/Lun/LunPressure.C
    kineticTheoryModels/granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C

    kineticTheoryModels/frictionalStressModel/frictionalStressModel/frictionalStressModel.C
    kineticTheoryModels/frictionalStressModel/frictionalStressModel/frictionalStressModelNew.C
    kineticTheoryModels/frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C

    derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.C
    derivedFvPatchFields/JohnsonJacksonParticleTheta/JohnsonJacksonParticleThetaFvPatchScalarField.C

)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE

    ${LIB_SRC}/meshTools/lnInclude
    ${LIB_SRC}/finiteVolume/lnInclude
    ${LIB_SRC}/sampling/lnInclude
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude
    ${LIB_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude
    ${LIB_SRC}/twoPhaseModels/compressibleTwoPhases/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/MomentumTransportModels/phaseCompressible/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    phaseSystem
    fluidThermophysicalTransportModel
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


