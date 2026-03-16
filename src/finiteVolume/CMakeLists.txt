lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_library(finiteVolume SHARED
    fvMesh/fvMeshGeometry.C
    fvMesh/fvMesh.C

    fvMesh/singleCellFvMesh/singleCellFvMesh.C

    fvMesh/zeroDimensionalFvMesh/zeroDimensionalFvMesh.C

    fvMesh/fvCellZone/fvCellZone.C
)


set(fvBoundaryMesh fvMesh/fvBoundaryMesh)
target_sources(finiteVolume PRIVATE
    ${fvBoundaryMesh}/fvBoundaryMesh.C
)

set(fvPatches fvMesh/fvPatches)
target_sources(finiteVolume PRIVATE
    ${fvPatches}/fvPatch/fvPatch.C
    ${fvPatches}/fvPatch/fvPatchNew.C

)

set(basicFvPatches ${fvPatches}/basic)
target_sources(finiteVolume PRIVATE
    ${basicFvPatches}/coupled/coupledFvPatch.C
)

set(constraintFvPatches ${fvPatches}/constraint)
target_sources(finiteVolume PRIVATE
    ${constraintFvPatches}/cyclic/cyclicFvPatch.C
    ${constraintFvPatches}/cyclicSlip/cyclicSlipFvPatch.C
    ${constraintFvPatches}/empty/emptyFvPatch.C
    ${constraintFvPatches}/processor/processorFvPatch.C
    ${constraintFvPatches}/processorCyclic/processorCyclicFvPatch.C
    ${constraintFvPatches}/symmetryPlane/symmetryPlaneFvPatch.C
    ${constraintFvPatches}/symmetry/symmetryFvPatch.C
    ${constraintFvPatches}/wedge/wedgeFvPatch.C
    ${constraintFvPatches}/internal/internalFvPatch.C
    ${constraintFvPatches}/nonConformal/nonConformalFvPatch.C
    ${constraintFvPatches}/nonConformalCoupled/nonConformalCoupledFvPatch.C
    ${constraintFvPatches}/nonConformalCyclic/nonConformalCyclicFvPatch.C
    ${constraintFvPatches}/nonConformalError/nonConformalErrorFvPatch.C
    ${constraintFvPatches}/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatch.C
)

set(derivedFvPatches ${fvPatches}/derived)
target_sources(finiteVolume PRIVATE
    ${derivedFvPatches}/wall/wallFvPatch.C
    ${derivedFvPatches}/mapped/mappedFvPatchBaseBase.C
    ${derivedFvPatches}/mapped/mappedFvPatchBase.C
    ${derivedFvPatches}/mapped/mappedFvPatch.C
    ${derivedFvPatches}/mapped/mappedWallFvPatch.C
    ${derivedFvPatches}/mapped/mappedExtrudedWallFvPatch.C
    ${derivedFvPatches}/mapped/mappedInternalFvPatch.C
    ${derivedFvPatches}/nonConformalMapped/nonConformalMappedFvPatchBase.C
    ${derivedFvPatches}/nonConformalMapped/nonConformalMappedWallFvPatch.C

)

set(meshWave algorithms/FvFaceCellWave)
target_sources(finiteVolume PRIVATE
    ${meshWave}/FvFaceCellWaveName.C
)

set(wallDist fvMesh/wallDist)
target_sources(finiteVolume PRIVATE
    ${wallDist}/fvPatchDistWave/fvPatchDistWave.C
    ${wallDist}/nearWallDist/nearWallDist.C
    ${wallDist}/wallDist/wallDist.C
    ${wallDist}/patchDistMethods/patchDistMethod/patchDistMethod.C
    ${wallDist}/patchDistMethods/meshWave/meshWavePatchDistMethod.C
    ${wallDist}/patchDistMethods/Poisson/PoissonPatchDistMethod.C
    ${wallDist}/patchDistMethods/advectionDiffusion/advectionDiffusionPatchDistMethod.C
)

set(fvMeshMapper fvMesh/fvMeshMapper)
target_sources(finiteVolume PRIVATE
    ${fvMeshMapper}/fvPatchMapper.C
    ${fvMeshMapper}/fvSurfaceMapper.C
)

set(extendedStencil fvMesh/extendedStencil)
set(cellToCell ${extendedStencil}/cellToCell)
target_sources(finiteVolume PRIVATE
    ${cellToCell}/extendedCellToCellStencil.C
    ${cellToCell}/extendedCentredCellToCellStencil.C
    ${cellToCell}/globalIndexStencils/cellToCellStencil.C
    ${cellToCell}/globalIndexStencils/CFCCellToCellStencil.C
    ${cellToCell}/globalIndexStencils/CPCCellToCellStencil.C
    ${cellToCell}/globalIndexStencils/CECCellToCellStencil.C
    ${cellToCell}/MeshObjects/centredCECCellToCellStencilObject.C
    ${cellToCell}/MeshObjects/centredCFCCellToCellStencilObject.C
    ${cellToCell}/MeshObjects/centredCPCCellToCellStencilObject.C
)

set(cellToFace ${extendedStencil}/cellToFace)
target_sources(finiteVolume PRIVATE
    ${cellToFace}/globalIndexStencils/cellToFaceStencil.C
    ${cellToFace}/globalIndexStencils/CFCCellToFaceStencil.C
    ${cellToFace}/globalIndexStencils/CECCellToFaceStencil.C
    ${cellToFace}/globalIndexStencils/CPCCellToFaceStencil.C
    ${cellToFace}/globalIndexStencils/FECCellToFaceStencil.C
    ${cellToFace}/extendedCellToFaceStencil.C
    ${cellToFace}/extendedCentredCellToFaceStencil.C
    ${cellToFace}/extendedUpwindCellToFaceStencil.C
    ${cellToFace}/MeshObjects/centredCECCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/centredCFCCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/centredCPCCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/centredFECCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/upwindCECCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/upwindCFCCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/upwindCPCCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/upwindFECCellToFaceStencilObject.C
    ${cellToFace}/MeshObjects/pureUpwindCFCCellToFaceStencilObject.C
)

set(faceToCell ${extendedStencil}/faceToCell)
target_sources(finiteVolume PRIVATE
    ${faceToCell}/globalIndexStencils/faceToCellStencil.C
    ${faceToCell}/globalIndexStencils/CFCFaceToCellStencil.C
    ${faceToCell}/extendedFaceToCellStencil.C
    ${faceToCell}/extendedCentredFaceToCellStencil.C
    ${faceToCell}/MeshObjects/centredCFCFaceToCellStencilObject.C
)

set(fvMeshTopoChangers fvMesh/fvMeshTopoChangers)
target_sources(finiteVolume PRIVATE
    ${fvMeshTopoChangers}/fvMeshTopoChanger/fvMeshTopoChanger.C
    ${fvMeshTopoChangers}/fvMeshTopoChanger/fvMeshTopoChangerNew.C
    ${fvMeshTopoChangers}/none/none_fvMeshTopoChanger.C
    ${fvMeshTopoChangers}/list/list_fvMeshTopoChanger.C
)

set(fvMeshDistributors fvMesh/fvMeshDistributors)
target_sources(finiteVolume PRIVATE
    ${fvMeshDistributors}/fvMeshDistributor/fvMeshDistributor.C
    ${fvMeshDistributors}/fvMeshDistributor/fvMeshDistributorNew.C
    ${fvMeshDistributors}/none/none_fvMeshDistributor.C
)

set(fvMeshMovers fvMesh/fvMeshMovers)
target_sources(finiteVolume PRIVATE
    ${fvMeshMovers}/fvMeshMover/fvMeshMover.C
    ${fvMeshMovers}/fvMeshMover/fvMeshMoverNew.C
    ${fvMeshMovers}/none/none_fvMeshMover.C
)

set(fvMeshStitchers fvMesh/fvMeshStitchers)
target_sources(finiteVolume PRIVATE
    ${fvMeshStitchers}/fvMeshStitcher/fvMeshStitcher.C
    ${fvMeshStitchers}/fvMeshStitcher/fvMeshStitcherNew.C
    ${fvMeshStitchers}/fvMeshStitcher/fvMeshStitcherTools.C
    ${fvMeshStitchers}/fvMeshStitcher/conformedFvPatchFields.C
    ${fvMeshStitchers}/fvMeshStitcher/conformedFvsPatchFields.C
    ${fvMeshStitchers}/stationary/stationary_fvMeshStitcher.C

    functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.C
    functionObjects/fvModelFunctionObject/fvModelFunctionObject.C
)

set(fvPatchFields fields/fvPatchFields)
target_sources(finiteVolume PRIVATE
    ${fvPatchFields}/fvPatchField/fvPatchFields.C

    ${fvPatchFields}/DimensionedFvPatchFields/DimensionedFvPatchFields.C
    ${fvPatchFields}/DimensionedFvPatchFieldFunctions/makeDimensionedFvPatchFieldFunctions.C
)

set(basicFvPatchFields ${fvPatchFields}/basic)
target_sources(finiteVolume PRIVATE
    ${basicFvPatchFields}/basicSymmetry/basicSymmetryFvPatchScalarField.C
    ${basicFvPatchFields}/calculated/calculatedFvPatchFields.C
    ${basicFvPatchFields}/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.C
    ${basicFvPatchFields}/coupled/coupledFvPatchFields.C
    ${basicFvPatchFields}/directionMixed/directionMixedFvPatchFields.C
    ${basicFvPatchFields}/fixedGradient/fixedGradientFvPatchFields.C
    ${basicFvPatchFields}/fixedValue/fixedValueFvPatchFields.C
    ${basicFvPatchFields}/mixed/mixedFvPatchFields.C
    ${basicFvPatchFields}/sliced/slicedFvPatchFields.C
    ${basicFvPatchFields}/transform/transformFvPatchFields.C
    ${basicFvPatchFields}/transform/transformFvPatchScalarField.C
    ${basicFvPatchFields}/zeroGradient/zeroGradientFvPatchFields.C
)

set(constraintFvPatchFields ${fvPatchFields}/constraint)
target_sources(finiteVolume PRIVATE
    ${constraintFvPatchFields}/cyclic/cyclicFvPatchFields.C
    ${constraintFvPatchFields}/cyclicSlip/cyclicSlipFvPatchFields.C
    ${constraintFvPatchFields}/empty/emptyFvPatchFields.C
    ${constraintFvPatchFields}/jumpCyclic/jumpCyclicFvPatchFields.C
    ${constraintFvPatchFields}/processor/processorFvPatchFields.C
    ${constraintFvPatchFields}/processor/processorFvPatchScalarField.C
    ${constraintFvPatchFields}/processorCyclic/processorCyclicFvPatchFields.C
    ${constraintFvPatchFields}/symmetryPlane/symmetryPlaneFvPatchFields.C
    ${constraintFvPatchFields}/symmetryPlane/symmetryPlaneFvPatchScalarField.C
    ${constraintFvPatchFields}/symmetry/symmetryFvPatchFields.C
    ${constraintFvPatchFields}/wedge/wedgeFvPatchFields.C
    ${constraintFvPatchFields}/wedge/wedgeFvPatchScalarField.C
    ${constraintFvPatchFields}/internal/internalFvPatchFields.C
    ${constraintFvPatchFields}/nonConformalCyclic/nonConformalCyclicFvPatchFields.C
    ${constraintFvPatchFields}/nonConformalError/nonConformalErrorFvPatchFields.C
    ${constraintFvPatchFields}/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatchFields.C

)

set(derivedFvPatchFields ${fvPatchFields}/derived)
target_sources(finiteVolume PRIVATE
    ${derivedFvPatchFields}/advective/advectiveFvPatchFields.C
    ${derivedFvPatchFields}/codedFixedValue/codedFixedValueFvPatchFields.C
    ${derivedFvPatchFields}/functionalFixedValue/functionalFixedValueFvPatchFields.C
    ${derivedFvPatchFields}/codedMixed/codedMixedFvPatchFields.C
    ${derivedFvPatchFields}/externalCoupledMixed/externalCoupledMixedFvPatchFields.C
    ${derivedFvPatchFields}/fanPressureJump/fanPressureJumpFvPatchScalarField.C
    ${derivedFvPatchFields}/fanPressure/fanPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/fixedFluxPressure/fixedFluxPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/fixedFluxExtrapolatedPressure/fixedFluxExtrapolatedPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/fixedInternalValue/fixedInternalValueFvPatchFields.C
    ${derivedFvPatchFields}/fixedJump/fixedJumpFvPatchFields.C
    ${derivedFvPatchFields}/fixedMean/fixedMeanFvPatchFields.C
    ${derivedFvPatchFields}/fixedNormalSlip/fixedNormalSlipFvPatchFields.C
    ${derivedFvPatchFields}/fixedPressureCompressibleDensity/fixedPressureCompressibleDensityFvPatchScalarField.C
    ${derivedFvPatchFields}/flowRateInletVelocity/flowRateInletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/flowRateInletVelocity/laminarBL/laminarBL.C
    ${derivedFvPatchFields}/flowRateInletVelocity/turbulentBL/turbulentBL.C
    ${derivedFvPatchFields}/flowRateOutletVelocity/flowRateOutletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/fluxCorrectedVelocity/fluxCorrectedVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/freestream/freestreamFvPatchFields.C
    ${derivedFvPatchFields}/freestreamVelocity/freestreamVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/freestreamPressure/freestreamPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/inletOutlet/inletOutletFvPatchFields.C
    ${derivedFvPatchFields}/zeroInletOutlet/zeroInletOutletFvPatchFields.C
    ${derivedFvPatchFields}/inletOutletTotalTemperature/inletOutletTotalTemperatureFvPatchScalarField.C
    ${derivedFvPatchFields}/mappedFlowRateVelocity/mappedFlowRateVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/mappedInternalValue/mappedInternalValueFvPatchFields.C
    ${derivedFvPatchFields}/mappedValue/mappedValueFvPatchFields.C
    ${derivedFvPatchFields}/mappedVelocityFlux/mappedVelocityFluxFvPatchField.C
    ${derivedFvPatchFields}/matchedFlowRateOutletVelocity/matchedFlowRateOutletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/noSlip/noSlipFvPatchVectorField.C
    ${derivedFvPatchFields}/movingWallVelocity/movingWallVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/movingWallSlipVelocity/movingWallSlipVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/movingMappedWallVelocity/movingMappedWallVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/outletInlet/outletInletFvPatchFields.C
    ${derivedFvPatchFields}/outletMappedUniformInlet/outletMappedUniformInletFvPatchFields.C
    ${derivedFvPatchFields}/fixedMeanOutletInlet/fixedMeanOutletInletFvPatchFields.C
    ${derivedFvPatchFields}/partialSlip/partialSlipFvPatchFields.C
    ${derivedFvPatchFields}/phaseHydrostaticPressure/phaseHydrostaticPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/pressureDirectedInletOutletVelocity/pressureDirectedInletOutletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/pressureDirectedInletVelocity/pressureDirectedInletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/pressureInletOutletParSlipVelocity/pressureInletOutletParSlipVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/pressureInletOutletVelocity/pressureInletOutletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/pressureInletUniformVelocity/pressureInletUniformVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/pressureInletVelocity/pressureInletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/pressureNormalInletOutletVelocity/pressureNormalInletOutletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/fixedNormalInletOutletVelocity/fixedNormalInletOutletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/rotatingPressureInletOutletVelocity/rotatingPressureInletOutletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/rotatingTotalPressure/rotatingTotalPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/rotatingWallVelocity/rotatingWallVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/slip/slipFvPatchFields.C
    ${derivedFvPatchFields}/supersonicFreestream/supersonicFreestreamFvPatchVectorField.C
    ${derivedFvPatchFields}/surfaceNormalFixedValue/surfaceNormalFixedValueFvPatchVectorField.C
    ${derivedFvPatchFields}/surfaceNormalUniformFixedValue/surfaceNormalUniformFixedValueFvPatchVectorField.C
    ${derivedFvPatchFields}/swirlFlowRateInletVelocity/swirlFlowRateInletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/swirlInletVelocity/swirlInletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/syringePressure/syringePressureFvPatchScalarField.C
    ${derivedFvPatchFields}/timeVaryingMappedFixedValue/timeVaryingMappedFvPatchFieldName.C
    ${derivedFvPatchFields}/timeVaryingMappedFixedValue/timeVaryingMappedFixedValueFvPatchFields.C
    ${derivedFvPatchFields}/dynamicPressure/dynamicPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/totalPressure/totalPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/entrainmentPressure/entrainmentPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/totalTemperature/totalTemperatureFvPatchScalarField.C
    ${derivedFvPatchFields}/translatingWallVelocity/translatingWallVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/turbulentInlet/turbulentInletFvPatchFields.C
    ${derivedFvPatchFields}/turbulentIntensityKineticEnergyInlet/turbulentIntensityKineticEnergyInletFvPatchScalarField.C
    ${derivedFvPatchFields}/uniformDensityHydrostaticPressure/uniformDensityHydrostaticPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/uniformFixedGradient/uniformFixedGradientFvPatchFields.C
    ${derivedFvPatchFields}/uniformFixedValue/uniformFixedValueFvPatchFields.C
    ${derivedFvPatchFields}/uniformInletOutlet/uniformInletOutletFvPatchFields.C
    ${derivedFvPatchFields}/uniformJump/uniformJumpFvPatchFields.C
    ${derivedFvPatchFields}/uniformTotalPressure/uniformTotalPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/variableHeightFlowRate/variableHeightFlowRateFvPatchField.C
    ${derivedFvPatchFields}/variableHeightFlowRateInletVelocity/variableHeightFlowRateInletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/outletPhaseMeanVelocity/outletPhaseMeanVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/waveTransmissive/waveTransmissiveFvPatchFields.C
    ${derivedFvPatchFields}/waveSurfacePressure/waveSurfacePressureFvPatchScalarField.C
    ${derivedFvPatchFields}/interstitialInletVelocity/interstitialInletVelocityFvPatchVectorField.C
    ${derivedFvPatchFields}/fixedProfile/fixedProfileFvPatchFields.C
    ${derivedFvPatchFields}/plenumPressure/plenumPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/interfaceCompression/interfaceCompressionFvPatchScalarField.C
    ${derivedFvPatchFields}/pressure/pressureFvPatchScalarField.C
    ${derivedFvPatchFields}/PrghPressure/prghPressureFvPatchScalarFields.C
    ${derivedFvPatchFields}/prghTotalHydrostaticPressure/prghTotalHydrostaticPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/fixedValueInletOutlet/fixedValueInletOutletFvPatchFields.C
    ${derivedFvPatchFields}/transonicEntrainmentPressure/transonicEntrainmentPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/prghCyclicPressure/prghCyclicPressureFvPatchScalarField.C
    ${derivedFvPatchFields}/zeroFixedValue/zeroFixedValueFvPatchFields.C
    ${derivedFvPatchFields}/meanInletOutlet/meanInletOutletFvPatchScalarField.C
)

set(fvsPatchFields fields/fvsPatchFields)
target_sources(finiteVolume PRIVATE
    ${fvsPatchFields}/fvsPatchField/fvsPatchFields.C
)

set(basicFvsPatchFields ${fvsPatchFields}/basic)
target_sources(finiteVolume PRIVATE
    ${basicFvsPatchFields}/calculated/calculatedFvsPatchFields.C
    ${basicFvsPatchFields}/coupled/coupledFvsPatchFields.C
    ${basicFvsPatchFields}/fixedValue/fixedValueFvsPatchFields.C
    ${basicFvsPatchFields}/sliced/slicedFvsPatchFields.C
)

set(constraintFvsPatchFields ${fvsPatchFields}/constraint)
target_sources(finiteVolume PRIVATE
    ${constraintFvsPatchFields}/cyclic/cyclicFvsPatchFields.C
    ${constraintFvsPatchFields}/cyclicSlip/cyclicSlipFvsPatchFields.C
    ${constraintFvsPatchFields}/empty/emptyFvsPatchFields.C
    ${constraintFvsPatchFields}/processor/processorFvsPatchFields.C
    ${constraintFvsPatchFields}/processorCyclic/processorCyclicFvsPatchFields.C
    ${constraintFvsPatchFields}/symmetryPlane/symmetryPlaneFvsPatchFields.C
    ${constraintFvsPatchFields}/symmetry/symmetryFvsPatchFields.C
    ${constraintFvsPatchFields}/wedge/wedgeFvsPatchFields.C
    ${constraintFvsPatchFields}/internal/internalFvsPatchFields.C
    ${constraintFvsPatchFields}/nonConformalCyclic/nonConformalCyclicFvsPatchFields.C
    ${constraintFvsPatchFields}/nonConformalError/nonConformalErrorFvsPatchFields.C
    ${constraintFvsPatchFields}/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvsPatchFields.C
)

set(derivedFvsPatchFields ${fvsPatchFields}/derived)
target_sources(finiteVolume PRIVATE
    ${derivedFvsPatchFields}/zeroFixedValue/zeroFixedValueFvsPatchFields.C
    ${derivedFvsPatchFields}/polyFaces/polyFacesFvsPatchLabelField.C
    ${derivedFvsPatchFields}/nonConformalPolyFaces/nonConformalPolyFacesFvsPatchLabelField.C
    ${derivedFvsPatchFields}/nonConformalMappedPolyFaces/nonConformalMappedPolyFacesFvsPatchLabelField.C

)

set(fvFieldSources fields/fvFieldSources)
target_sources(finiteVolume PRIVATE
    ${fvFieldSources}/fvFieldSource/fvFieldSources.C
)

set(derivedFvFieldSources ${fvFieldSources}/derived)
target_sources(finiteVolume PRIVATE
    ${derivedFvFieldSources}/internal/internalFvFieldSources.C
    ${derivedFvFieldSources}/uniformFixedValue/uniformFixedValueFvFieldSources.C
    ${derivedFvFieldSources}/uniformInletOutlet/uniformInletOutletFvFieldSources.C
    ${derivedFvFieldSources}/turbulentIntensityKineticEnergy/turbulentIntensityKineticEnergyFvScalarFieldSource.C
    ${derivedFvFieldSources}/NaN/NaNFvFieldSources.C
)

target_sources(finiteVolume PRIVATE
    fields/DimensionedFields/DimensionedFieldFunctions/makeDimensionedVolFieldFunctions.C
    fields/DimensionedFields/DimensionedFieldFunctions/makeDimensionedSurfaceFieldFunctions.C

    fields/volFields/volFields.C
    fields/surfaceFields/surfaceFields.C

    fvMatrices/fvMatrices.C
    fvMatrices/fvScalarMatrix/fvScalarMatrix.C
    fvMatrices/solvers/MULES/MULES.C
    fvMatrices/solvers/GAMGSymSolver/GAMGAgglomerations/faceAreaPairGAMGAgglomeration/faceAreaPairGAMGAgglomeration.C

)

set(interpolation interpolation/interpolation)
target_sources(finiteVolume PRIVATE
    ${interpolation}/interpolation/interpolations.C
    ${interpolation}/cell/cell_interpolations.C
    ${interpolation}/cellPatchConstrained/cellPatchConstraineds.C
    ${interpolation}/volPointInterpolation/volPointInterpolation_interpolations.C
    ${interpolation}/cellPoint/cellPoint_interpolations.C
    ${interpolation}/cellPointFace/cellPointFaces.C
    ${interpolation}/cellPointWallModified/cellPointWallModifieds.C
    ${interpolation}/pointMVC/pointMVCWeight.C
    ${interpolation}/pointMVC/pointMVCs.C
)

set(volPointInterpolation interpolation/volPointInterpolation)
target_sources(finiteVolume PRIVATE
    ${volPointInterpolation}/volPointInterpolation.C
    ${volPointInterpolation}/pointConstraints.C
)

set(surfaceInterpolation interpolation/surfaceInterpolation)
target_sources(finiteVolume PRIVATE
    ${surfaceInterpolation}/surfaceInterpolation/surfaceInterpolation.C
    ${surfaceInterpolation}/surfaceInterpolationScheme/surfaceInterpolationSchemes.C

    ${surfaceInterpolation}/blendedSchemeBase/blendedSchemeBaseName.C
)

set(schemes ${surfaceInterpolation}/schemes)
target_sources(finiteVolume PRIVATE
    ${schemes}/linear/linear.C
    ${schemes}/pointLinear/pointLinear.C
    ${schemes}/midPoint/midPoint.C
    ${schemes}/downwind/downwind.C
    ${schemes}/weighted/weighted.C
    ${schemes}/cubic/cubic.C
    ${schemes}/skewCorrected/skewCorrectionVectors.C
    ${schemes}/skewCorrected/skewCorrected.C
    ${schemes}/outletStabilised/outletStabilised.C
    ${schemes}/phaseStabilised/phaseStabilised.C
    ${schemes}/reverseLinear/reverseLinear.C
    ${schemes}/clippedLinear/clippedLinear.C
    ${schemes}/harmonic/harmonic.C
    ${schemes}/fixedBlended/fixedBlended.C
    ${schemes}/localBlended/localBlended.C
    ${schemes}/limiterBlended/limiterBlended.C
    ${schemes}/CoBlended/CoBlended.C
    ${schemes}/cellCoBlended/cellCoBlended.C
    ${schemes}/deferred/deferred.C
    ${schemes}/localMax/localMax.C
    ${schemes}/localMin/localMin.C

    ${schemes}/linearFit/linearFit.C
    ${schemes}/biLinearFit/biLinearFit.C
    ${schemes}/quadraticLinearFit/quadraticLinearFit.C
    ${schemes}/quadraticFit/quadraticFit.C

    ${schemes}/quadraticLinearUpwindFit/quadraticLinearUpwindFit.C
    ${schemes}/quadraticUpwindFit/quadraticUpwindFit.C
    ${schemes}/cubicUpwindFit/cubicUpwindFit.C
    # ${schemes}/quadraticLinearPureUpwindFit/quadraticLinearPureUpwindFit.C
    ${schemes}/linearPureUpwindFit/linearPureUpwindFit.C
    ${schemes}/linearUpwind/linearUpwind.C
    ${schemes}/linearUpwind/linearUpwindV.C
    ${schemes}/LUST/LUST.C
)

set(limitedSchemes ${surfaceInterpolation}/limitedSchemes)
target_sources(finiteVolume PRIVATE
    ${limitedSchemes}/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationSchemes.C
    ${limitedSchemes}/upwind/upwind.C
    ${limitedSchemes}/blended/blended.C
    ${limitedSchemes}/Gamma/Gamma.C
    ${limitedSchemes}/SFCD/SFCD.C
    ${limitedSchemes}/Minmod/Minmod.C
    ${limitedSchemes}/vanLeer/vanLeer.C
    ${limitedSchemes}/vanAlbada/vanAlbada.C
    ${limitedSchemes}/OSPRE/OSPRE.C
    ${limitedSchemes}/limitedLinear/limitedLinear.C
    ${limitedSchemes}/limitedCubic/limitedCubic.C
    ${limitedSchemes}/SuperBee/SuperBee.C
    ${limitedSchemes}/QUICK/QUICK.C
    ${limitedSchemes}/MUSCL/MUSCL.C
    ${limitedSchemes}/UMIST/UMIST.C
    ${limitedSchemes}/Phi/Phi.C
    ${limitedSchemes}/filteredLinear/filteredLinear.C
    ${limitedSchemes}/filteredLinear2/filteredLinear2.C
    ${limitedSchemes}/filteredLinear3/filteredLinear3.C
    ${limitedSchemes}/limitWith/limitWith.C
)

set(multivariateSchemes ${surfaceInterpolation}/multivariateSchemes)
target_sources(finiteVolume PRIVATE
    ${multivariateSchemes}/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationSchemes.C
    ${multivariateSchemes}/multivariateSelectionScheme/multivariateSelectionSchemes.C
    ${multivariateSchemes}/multivariateIndependentScheme/multivariateIndependentSchemes.C
    ${multivariateSchemes}/upwind/multivariateUpwind.C
    ${multivariateSchemes}/Gamma/multivariateGamma.C
    ${multivariateSchemes}/vanLeer/multivariateVanLeer.C
    ${multivariateSchemes}/Minmod/multivariateMinmod.C
    ${multivariateSchemes}/SuperBee/multivariateSuperBee.C
    ${multivariateSchemes}/MUSCL/multivariateMUSCL.C
    ${multivariateSchemes}/limitedLinear/multivariateLimitedLinear.C
    ${multivariateSchemes}/limitedCubic/multivariateLimitedCubic.C
)

target_sources(finiteVolume PRIVATE
    finiteVolume/fv/fv.C
    finiteVolume/fvSchemes/fvSchemes.C
)

set(ddtSchemes finiteVolume/ddtSchemes)
target_sources(finiteVolume PRIVATE
    ${ddtSchemes}/ddtScheme/ddtSchemes.C
    ${ddtSchemes}/steadyStateDdtScheme/steadyStateDdtSchemes.C
    ${ddtSchemes}/EulerDdtScheme/EulerDdtSchemes.C
    ${ddtSchemes}/CoEulerDdtScheme/CoEulerDdtSchemes.C
    ${ddtSchemes}/SLTSDdtScheme/SLTSDdtSchemes.C
    ${ddtSchemes}/localEulerDdtScheme/localEulerDdt.C
    ${ddtSchemes}/localEulerDdtScheme/localEulerDdtSchemes.C
    ${ddtSchemes}/backwardDdtScheme/backwardDdtSchemes.C
    ${ddtSchemes}/CrankNicolsonDdtScheme/CrankNicolsonDdtSchemes.C
    ${ddtSchemes}/boundedDdtScheme/boundedDdtSchemes.C
)

set(d2dt2Schemes finiteVolume/d2dt2Schemes)
target_sources(finiteVolume PRIVATE
    ${d2dt2Schemes}/d2dt2Scheme/d2dt2Schemes.C
    ${d2dt2Schemes}/steadyStateD2dt2Scheme/steadyStateD2dt2Schemes.C
    ${d2dt2Schemes}/EulerD2dt2Scheme/EulerD2dt2Schemes.C
)

set(divSchemes finiteVolume/divSchemes)
target_sources(finiteVolume PRIVATE
    ${divSchemes}/divScheme/divSchemes.C
    ${divSchemes}/gaussDivScheme/gaussDivSchemes.C
)

set(gradSchemes finiteVolume/gradSchemes)
target_sources(finiteVolume PRIVATE
    ${gradSchemes}/gradScheme/gradSchemes.C
    ${gradSchemes}/gaussGrad/gaussGrads.C

    ${gradSchemes}/leastSquaresGrad/leastSquaresVectors.C
    ${gradSchemes}/leastSquaresGrad/leastSquaresGrads.C
    ${gradSchemes}/LeastSquaresGrad/LeastSquaresGrads.C
    ${gradSchemes}/fourthGrad/fourthGrads.C
)

set(limitedGradSchemes ${gradSchemes}/limitedGradSchemes)
target_sources(finiteVolume PRIVATE
    ${limitedGradSchemes}/faceLimitedGrad/faceLimitedGrads.C
    ${limitedGradSchemes}/cellLimitedGrad/cellLimitedGrads.C
    ${limitedGradSchemes}/faceMDLimitedGrad/faceMDLimitedGrads.C
    ${limitedGradSchemes}/cellMDLimitedGrad/cellMDLimitedGrads.C
)

set(snGradSchemes finiteVolume/snGradSchemes)
target_sources(finiteVolume PRIVATE
    ${snGradSchemes}/snGradScheme/snGradSchemes.C
    ${snGradSchemes}/correctedSnGrad/correctedSnGrads.C
    ${snGradSchemes}/faceCorrectedSnGrad/faceCorrectedSnGrads.C
    ${snGradSchemes}/limitedSnGrad/limitedSnGrads.C
    ${snGradSchemes}/uncorrectedSnGrad/uncorrectedSnGrads.C
    ${snGradSchemes}/orthogonalSnGrad/orthogonalSnGrads.C
    ${snGradSchemes}/quadraticFitSnGrad/quadraticFitSnGrads.C
    ${snGradSchemes}/linearFitSnGrad/linearFitSnGrads.C
    ${snGradSchemes}/phaseStabilisedSnGrad/phaseStabilisedSnGrads.C
)

set(convectionSchemes finiteVolume/convectionSchemes)
target_sources(finiteVolume PRIVATE
    ${convectionSchemes}/convectionScheme/convectionSchemes.C
    ${convectionSchemes}/gaussConvectionScheme/gaussConvectionSchemes.C
    ${convectionSchemes}/multivariateGaussConvectionScheme/multivariateGaussConvectionSchemes.C
    ${convectionSchemes}/boundedConvectionScheme/boundedConvectionSchemes.C
)

set(laplacianSchemes finiteVolume/laplacianSchemes)
target_sources(finiteVolume PRIVATE
    ${laplacianSchemes}/laplacianScheme/laplacianSchemes.C
    ${laplacianSchemes}/gaussLaplacianScheme/gaussLaplacianSchemes.C
)

target_sources(finiteVolume PRIVATE
    finiteVolume/fvc/fvcMeshPhi.C
    finiteVolume/fvc/fvcSmooth/fvcSmooth.C
    finiteVolume/fvc/fvcReconstructMag.C
)

set(general cfdTools/general)
target_sources(finiteVolume PRIVATE
    ${general}/findRefCell/findRefCell.C
    ${general}/constrainHbyA/constrainHbyA.C
    ${general}/adjustPhi/adjustPhi.C
    ${general}/bound/bound.C
    ${general}/correctPhi/fvCorrectPhi.C
    ${general}/correctPhi/correctUphiBCs.C
    ${general}/pressureReference/pressureReference.C
    ${general}/levelSet/levelSet.C
    ${general}/surfaceToVolVelocity/surfaceToVolVelocity.C
    ${general}/buoyancy/buoyancy.C
)

set(solutionControl ${general}/solutionControl)
target_sources(finiteVolume PRIVATE
    ${solutionControl}/solutionControl/singleRegionSolutionControl/singleRegionSolutionControl.C
    ${solutionControl}/solutionControl/multiRegionSolutionControl/multiRegionSolutionControl.C
    ${solutionControl}/solutionControl/nonOrthogonalSolutionControl/nonOrthogonalSolutionControl.C
    ${solutionControl}/solutionControl/fluidSolutionControl/fluidSolutionControl.C
    ${solutionControl}/simpleControl/simpleControl.C
    ${solutionControl}/pisoControl/pisoControl.C
)

set(pimpleControl ${solutionControl}/pimpleControl)
target_sources(finiteVolume PRIVATE
    ${pimpleControl}/pimpleControl/pimpleControl.C
    ${pimpleControl}/pimpleNoLoopControl/pimpleNoLoopControl.C
    ${pimpleControl}/pimpleLoop/pimpleLoop.C
    ${pimpleControl}/pimpleSingleRegionControl/pimpleSingleRegionControl.C
    ${pimpleControl}/pimpleMultiRegionControl/pimpleMultiRegionControl.C
)

set(convergenceControl ${solutionControl}/convergenceControl)
target_sources(finiteVolume PRIVATE
    ${convergenceControl}/convergenceControl/convergenceControl.C
    ${convergenceControl}/singleRegionConvergenceControl/singleRegionConvergenceControl.C
    ${convergenceControl}/correctorConvergenceControl/correctorConvergenceControl.C
    ${convergenceControl}/singleRegionCorrectorConvergenceControl/singleRegionCorrectorConvergenceControl.C
)

set(porosity ${general}/porosityModel)
target_sources(finiteVolume PRIVATE
    ${porosity}/porosityModel/porosityModel.C
    ${porosity}/porosityModel/porosityModelNew.C
    ${porosity}/porosityModel/porosityModelList.C
    ${porosity}/porosityModel/IOporosityModelList.C
    ${porosity}/DarcyForchheimer/DarcyForchheimer.C
    ${porosity}/fixedCoeff/fixedCoeff.C
    ${porosity}/powerLaw/powerLaw.C
    ${porosity}/solidification/solidification.C
)

set(MRF ${general}/MRF)
target_sources(finiteVolume PRIVATE
    ${MRF}/MRFZone.C
    ${MRF}/MRFZoneList.C
    ${MRF}/IOMRFZoneList.C
    ${MRF}/derivedFvPatchFields/MRFPatchField/MRFPatchField.C
    ${MRF}/derivedFvPatchFields/MRFnoSlip/MRFnoSlipFvPatchVectorField.C
    ${MRF}/derivedFvPatchFields/MRFslip/MRFslipFvPatchVectorField.C
    ${MRF}/derivedFvPatchFields/MRFFreestreamVelocity/MRFFreestreamVelocityFvPatchVectorField.C

)

set(fvModels ${general}/fvModels)
target_sources(finiteVolume PRIVATE
    ${fvModels}/fvModel.C
    ${fvModels}/fvModels.C
    ${fvModels}/includeFvModelEntry/includeFvModelEntry.C
)

set(fvConstraints ${general}/fvConstraints)
target_sources(finiteVolume PRIVATE
    ${fvConstraints}/fvConstraint.C
    ${fvConstraints}/fvConstraints.C
    ${fvConstraints}/includeFvConstraintEntry/includeFvConstraintEntry.C
)

target_sources(finiteVolume PRIVATE
    ${general}/fvSource/fvSource.C
    ${general}/fvSource/fvSpecificSource.C
    ${general}/fvSource/fvTotalSource.C
)

target_sources(finiteVolume PRIVATE
    solver/solver.C
    solver/solverNew.C
)

target_sources(finiteVolume PRIVATE

    fvMeshToFvMesh/fvMeshToFvMesh.C
)


target_include_directories(finiteVolume PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_link_libraries(finiteVolume PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
    triSurface
    meshTools
    Pstream_best
)


target_compile_definitions(finiteVolume PRIVATE
    LIB_NAME=libfiniteVolume.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};finiteVolume" PARENT_SCOPE)


