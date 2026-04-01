lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name XiFluidSolver)

add_library(${target_name} SHARED
    uRhoMulticomponentThermo/mixtures/uMixture/uMixture.C
    uRhoMulticomponentThermo/mixtures/uHomogeneousMixture/uHomogeneousMixture.C
    uRhoMulticomponentThermo/mixtures/uInhomogeneousMixture/uInhomogeneousMixture.C
    uRhoMulticomponentThermo/mixtures/uInhomogeneousEGRMixture/uInhomogeneousEGRMixture.C
    uRhoMulticomponentThermo/mixtures/uMulticomponentMixture/uMulticomponentMixture.C

    bRhoMulticomponentThermo/mixtures/bMixture/bMixture.C
    bRhoMulticomponentThermo/mixtures/bHomogeneousMixture/bHomogeneousMixture.C
    bRhoMulticomponentThermo/mixtures/bInhomogeneousMixture/bInhomogeneousMixture.C
    bRhoMulticomponentThermo/mixtures/bMulticomponentMixture/bMulticomponentMixture.C

    ubMixtureMaps/ubMixtureMap/ubMixtureMap.C
    ubMixtureMaps/uHomogeneousbHomogeneous/uHomogeneousbHomogeneous.C
    ubMixtureMaps/uInhomogeneousbInhomogeneous/uInhomogeneousbInhomogeneous.C
    ubMixtureMaps/uInhomogeneousEGRbInhomogeneous/uInhomogeneousEGRbInhomogeneous.C
    ubMixtureMaps/uMulticomponentbInhomogeneous/uMulticomponentbInhomogeneous.C
    ubMixtureMaps/uInhomogeneousbMulticomponent/uInhomogeneousbMulticomponent.C

    uRhoMulticomponentThermo/uRhoMulticomponentThermo.C
    uRhoMulticomponentThermo/uRhoMulticomponentThermos.C

    bRhoMulticomponentThermo/bRhoMulticomponentThermo.C
    bRhoMulticomponentThermo/bRhoMulticomponentThermos.C

    ubRhoThermo/ubRhoThermo.C
    ubRhoThermo/derivedFvPatchFields/ubCoupledTemperature/ubCoupledTemperatureFvPatchScalarField.C

    ubMomentumTransportModel/ubMomentumTransportModel.C

    XiFluid.C
    thermophysicalPredictor.C

    laminarFlameSpeed/laminarFlameSpeed/laminarFlameSpeed.C
    laminarFlameSpeed/laminarFlameSpeed/laminarFlameSpeedNew.C
    laminarFlameSpeed/constant/constant_laminarFlameSpeed.C
    laminarFlameSpeed/Gulder/Gulder_laminarFlameSpeed.C
    laminarFlameSpeed/RaviPetersen/RaviPetersen.C

    SuModels/SuModel/SuModel.C
    SuModels/SuModel/SuModelNew.C
    SuModels/uniformConstant/uniformConstant_SuModel.C
    SuModels/unstrained/unstrained.C
    SuModels/linearEquilibrium/linearEquilibrium_SuModel.C
    SuModels/transport/transport_SuModel.C

    XiModels/XiModel/XiModel.C
    XiModels/XiModel/XiModelNew.C
    XiModels/uniformConstant/uniformConstant.C
    XiModels/equilibrium/equilibrium.C
    XiModels/transport/transport.C

    XiModels/XiEqModels/XiEqModel/XiEqModel.C
    XiModels/XiEqModels/XiEqModel/XiEqModelNew.C
    XiModels/XiEqModels/constant/constant.C
    XiModels/XiEqModels/Gulder/Gulder.C
    XiModels/XiEqModels/instability/instability.C

    XiModels/XiProfiles/XiProfile/XiProfile.C
    XiModels/XiProfiles/XiProfile/XiProfileNew.C
    XiModels/XiProfiles/uniform/uniform.C
    XiModels/XiProfiles/linear/linear.C
    XiModels/XiProfiles/cubic/cubic.C

    XiModels/XiGModels/XiGModel/XiGModel.C
    XiModels/XiGModels/XiGModel/XiGModelNew.C
    XiModels/XiGModels/KTS/KTS.C
    XiModels/XiGModels/instability/instability_XiGModel.C

    fvModels/ignition/kernelShapes/kernelShape/kernelShape.C
    fvModels/ignition/kernelShapes/kernelShape/kernelShapeNew.C
    fvModels/ignition/kernelShapes/planar/planarKernelShape.C
    fvModels/ignition/kernelShapes/cylindrical/cylindricalKernelShape.C
    fvModels/ignition/kernelShapes/spherical/sphericalKernelShape.C

    fvModels/ignition/bXiIgnition/bXiIgnition.C
    fvModels/ignition/bXiTimedIgnition/bXiTimedIgnition.C
    fvModels/ignition/bXiIgnited/bXiIgnited.C
    fvModels/ignition/constantbXiIgnition/constantbXiIgnition.C
    fvModels/ignition/GaussianbXiIgnition/GaussianbXiIgnition.C
    fvModels/ignition/bXiKernelCorr/bXiKernelCorr.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluidThermo/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/coupledThermophysicalTransportModels/lnInclude
    ${LIB_SRC}/thermophysicalModels/specie/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/thermophysicalModels/multicomponentThermo/lnInclude

)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    basicFluidSolver 
    isothermalFluidSolver

    fluidThermophysicalTransportModel
    momentumTransportModels 
    compressibleMomentumTransportModels 
    phaseCompressibleMomentumTransportModels 
    phaseFluidThermophysicalTransportModels 
    coupledThermophysicalTransportModels 
    multicomponentThermophysicalModels
    reactionModels
    specie 
    finiteVolume 
    fvModels 
    fvConstraints 
    sampling 
    meshTools
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


