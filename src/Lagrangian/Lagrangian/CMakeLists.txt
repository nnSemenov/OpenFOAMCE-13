lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(Lagrangian SHARED
    fields/makeLagrangianFieldFunctions.C

    fields/labelIODynamicField.C
    fields/barycentricIODynamicField.C

    stateModel/stateModel.C
)

target_sources(Lagrangian PRIVATE
    # State
    LagrangianState/LagrangianState.C

    # Meshes
    LagrangianMesh/LagrangianMesh.C
    LagrangianMesh/LagrangianSubMesh/LagrangianSubMesh.C
    LagrangianMesh/LagrangianBoundaryMesh/LagrangianBoundaryMesh.C

    # Equations
    LagrangianEqn/LagrangianVectorSp.C
)

# Patches
set(LagrangianPatches LagrangianMesh/LagrangianPatches)
target_sources(Lagrangian PRIVATE
    ${LagrangianPatches}/LagrangianPatch/LagrangianPatch.C
)

set(basicLagrangianPatches ${LagrangianPatches}/basic)
set(constraintLagrangianPatches ${LagrangianPatches}/constraint)
target_sources(Lagrangian PRIVATE
    ${constraintLagrangianPatches}/empty/emptyLagrangianPatch.C
    ${constraintLagrangianPatches}/wedge/wedgeLagrangianPatch.C
    ${constraintLagrangianPatches}/symmetry/symmetryLagrangianPatch.C
    ${constraintLagrangianPatches}/symmetryPlane/symmetryPlaneLagrangianPatch.C
    ${constraintLagrangianPatches}/cyclic/cyclicLagrangianPatch.C
    ${constraintLagrangianPatches}/processor/processorLagrangianPatch.C
    ${constraintLagrangianPatches}/processorCyclic/processorCyclicLagrangianPatch.C
    ${constraintLagrangianPatches}/nonConformalCyclic/nonConformalCyclicLagrangianPatch.C
    ${constraintLagrangianPatches}/nonConformalError/nonConformalErrorLagrangianPatch.C
    ${constraintLagrangianPatches}/nonConformalProcessorCyclic/nonConformalProcessorCyclicLagrangianPatch.C
    ${constraintLagrangianPatches}/internal/internalLagrangianPatch.C
)

set(derivedLagrangianPatches ${LagrangianPatches}/derived)
target_sources(Lagrangian PRIVATE
    ${derivedLagrangianPatches}/wall/wallLagrangianPatch.C
)

# Patch Fields
set(LagrangianPatchFields fields/LagrangianPatchFields)
target_sources(Lagrangian PRIVATE
    ${LagrangianPatchFields}/LagrangianPatchField/LagrangianPatchFields.C
)

set(basicLagrangianPatchFields ${LagrangianPatchFields}/basic)
target_sources(Lagrangian PRIVATE
    ${basicLagrangianPatchFields}/calculated/calculatedLagrangianPatchFields.C
)

set(constraintLagrangianPatchFields ${LagrangianPatchFields}/constraint)
target_sources(Lagrangian PRIVATE
    ${constraintLagrangianPatchFields}/empty/emptyLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/wedge/wedgeLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/symmetry/symmetryLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/symmetryPlane/symmetryPlaneLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/cyclic/cyclicLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/processor/processorLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/processorCyclic/processorCyclicLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/nonConformalCyclic/nonConformalCyclicLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/nonConformalError/nonConformalErrorLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/nonConformalProcessorCyclic/nonConformalProcessorCyclicLagrangianPatchFields.C
    ${constraintLagrangianPatchFields}/internal/internalLagrangianPatchFields.C
)

set(derivedLagrangianPatchFields ${LagrangianPatchFields}/derived)
# Field Sources
set(LagrangianFieldSources fields/LagrangianFieldSources)
target_sources(Lagrangian PRIVATE
    ${LagrangianFieldSources}/LagrangianFieldSource/LagrangianFieldSourceBase.C
    ${LagrangianFieldSources}/LagrangianFieldSource/LagrangianFieldSources.C
)

set(derivedLagrangianFieldSources ${LagrangianFieldSources}/derived)
target_sources(Lagrangian PRIVATE
    ${derivedLagrangianFieldSources}/internal/internalLagrangianFieldSources.C
    ${derivedLagrangianFieldSources}/uniformFixedValue/uniformFixedValueLagrangianFieldSources.C
    ${derivedLagrangianFieldSources}/zero/zeroLagrangianFieldSources.C
    ${derivedLagrangianFieldSources}/max/maxLagrangianFieldSources.C
    ${derivedLagrangianFieldSources}/NaN/NaNLagrangianFieldSources.C
    ${derivedLagrangianFieldSources}/distribution/distributionLagrangianScalarFieldSource.C
    ${derivedLagrangianFieldSources}/noneState/noneStateLagrangianLabelFieldSource.C
    ${derivedLagrangianFieldSources}/position/positionLagrangianVectorFieldSource.C
)

target_sources(Lagrangian PRIVATE
    # Fields
    fields/LagrangianFields/LagrangianFields.C
    fields/LagrangianSubFields/LagrangianSubFields.C
    # Schemes
    Lagrangian/LagrangianSchemes/LagrangianSchemes.C
    Lagrangian/LagrangianDdtSchemes/LagrangianDdtScheme/LagrangianDdtSchemes.C
    Lagrangian/LagrangianDdtSchemes/Euler/Euler_LagrangianDdtSchemes.C
    Lagrangian/LagrangianDdtSchemes/CrankNicolson/CrankNicolson_LagrangianDdtSchemes.C
    Lagrangian/LagrangianSpSchemes/LagrangianSpScheme/LagrangianSpSchemes.C
    Lagrangian/LagrangianSpSchemes/implicit/implicit_LagrangianSpSchemes.C
    Lagrangian/LagrangianSpSchemes/Explicit/Explicit_LagrangianSpSchemes.C
    Lagrangian/LagrangianSpSchemes/upwind/upwind_LagrangianSpSchemes.C
    Lagrangian/LagrangianAccumulationSchemes/LagrangianAccumulationScheme/LagrangianAccumulationSchemes.C
    Lagrangian/LagrangianAccumulationSchemes/cell/cell_LagrangianAccumulationSchemes.C
    Lagrangian/LagrangianAccumulationSchemes/cellPoint/cellPoint_LagrangianAccumulationSchemes.C
    Lagrangian/LagrangianAccumulationSchemes/cellPoint/cellPointLagrangianAddressor.C
    Lagrangian/LagrangianAccumulationSchemes/cellPoint/cellPointLagrangianAccumulator.C

    # Solution
    Lagrangian/LagrangianSolution/LagrangianSolution.C

    # Averages
    LagrangianAverage/LagrangianAverage/LagrangianAverages.C
    LagrangianAverage/cell/cell_LagrangianAverages.C
    LagrangianAverage/cellPoint/cellPoint_LagrangianAverages.C

    # Models
    LagrangianModels/LagrangianModel/LagrangianModel.C
    LagrangianModels/LagrangianModels/LagrangianModels.C
    LagrangianModels/LagrangianSource/LagrangianSource.C
    LagrangianModels/LagrangianInjection/LagrangianInjection.C

    # Function objects
    LagrangianMeshFunctionObject/LagrangianMeshFunctionObject.C
)

target_include_directories(Lagrangian PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(Lagrangian PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    meshTools
    tracking
)

target_compile_definitions(Lagrangian PRIVATE
    LIB_NAME=libLagrangian.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};Lagrangian" PARENT_SCOPE)


