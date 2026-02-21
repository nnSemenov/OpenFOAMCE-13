
add_library(OpenFOAM SHARED
    global/argList/argList.C
    global/clock/clock.C
    global/etcFiles/etcFiles.C
)

set(processed_file ${CMAKE_CURRENT_BINARY_DIR}/global/global.C)
process_Cver(global/global.Cver ${processed_file})
target_sources(OpenFOAM PRIVATE
    ${processed_file}
)

set(fileOps global/fileOperations)
target_sources(OpenFOAM PRIVATE
    ${fileOps}/fileOperation/fileOperation.C
    ${fileOps}/fileOperationInitialise/fileOperationInitialise.C
    ${fileOps}/uncollatedFileOperation/uncollatedFileOperation.C
    ${fileOps}/masterUncollatedFileOperation/masterUncollatedFileOperation.C
    ${fileOps}/collatedFileOperation/collatedFileOperation.C
    ${fileOps}/collatedFileOperation/hostCollatedFileOperation.C
    ${fileOps}/collatedFileOperation/threadedCollatedOFstream.C
    ${fileOps}/collatedFileOperation/OFstreamCollator.C
)

set(bools primitives/bools)
target_sources(OpenFOAM PRIVATE
    ${bools}/bool/bool.C
    ${bools}/bool/boolIO.C
    ${bools}/Switch/Switch.C
    ${bools}/Switch/SwitchIO.C
    ${bools}/lists/boolList.C
)
set(chars primitives/chars)
target_sources(OpenFOAM PRIVATE
    ${chars}/char/charIO.C
    ${chars}/wchar/wcharIO.C
)
set(ints primitives/ints)
target_sources(OpenFOAM PRIVATE
    ${ints}/uint32/uint32.C
    ${ints}/uint32/uint32IO.C
    ${ints}/uint64/uint64.C
    ${ints}/uint64/uint64IO.C
    ${ints}/uint/uintIO.C
    ${ints}/int32/int32.C
    ${ints}/int32/int32IO.C
    ${ints}/int64/int64.C
    ${ints}/int64/int64IO.C
    ${ints}/int/intIO.C
    ${ints}/label/label.C
    ${ints}/uLabel/uLabel.C
    ${ints}/lists/labelIOList.C
    ${ints}/lists/labelListIOList.C
)
target_sources(OpenFOAM PRIVATE
    primitives/Scalar/longDoubleScalar/longDoubleScalar.C
    primitives/Scalar/doubleScalar/doubleScalar.C
    primitives/Scalar/floatScalar/floatScalar.C
    primitives/Scalar/scalar/scalar.C
    primitives/Scalar/scalar/incGamma.C
    primitives/Scalar/scalar/invIncGamma.C
    primitives/Scalar/lists/scalarList.C
    primitives/Scalar/lists/scalarIOList.C
    primitives/Scalar/lists/scalarListIOList.C

    primitives/DiagTensor/diagTensor/diagTensor.C
    primitives/SphericalTensor/sphericalTensor/sphericalTensor.C
    primitives/SphericalTensor/labelSphericalTensor/labelSphericalTensor.C
    primitives/SymmTensor/labelSymmTensor/labelSymmTensor.C
    primitives/SymmTensor/symmTensor/symmTensor.C
    primitives/Tensor/labelTensor/labelTensor.C
    primitives/Tensor/tensor/tensor.C
    primitives/Tensor/lists/sphericalTensorList.C
    primitives/Tensor/lists/symmTensorList.C
    primitives/Tensor/lists/tensorList.C
    primitives/direction/directionIO.C

    primitives/Vector/complexVector/complexVector.C
)

if (NOT ${WM_PRECISION_OPTION} STREQUAL "SP")
    target_sources(OpenFOAM PRIVATE
        primitives/Vector/floatVector/floatVector.C
        primitives/Tensor/floatTensor/floatTensor.C
    )
endif ()

target_sources(OpenFOAM PRIVATE
    primitives/Vector/labelVector/labelVector.C
    primitives/Vector/vector/vector.C
    primitives/Vector/lists/vectorList.C
    primitives/Vector/lists/vectorIOList.C
    primitives/Vector/lists/vectorListIOList.C

    primitives/Tensor2D/tensor2D/tensor2D.C
    primitives/SphericalTensor2D/sphericalTensor2D/sphericalTensor2D.C
    primitives/SymmTensor2D/symmTensor2D/symmTensor2D.C
    primitives/Vector2D/vector2D/vector2D.C

    primitives/complex/complex.C
    primitives/globalIndexAndTransform/globalIndexAndTransform.C
    primitives/transform/transformer/transformer.C
    primitives/transform/transformer/transformerIOList.C
    primitives/quaternion/quaternion.C
    primitives/septernion/septernion.C
    primitives/triad/triad.C

    primitives/omega/omega.C

    primitives/Tuple2/Tuple2s.C
)


# Run-time selectable functions
target_sources(OpenFOAM PRIVATE
    primitives/functions/Function1/Function1/Function1UnitConversions.C
    primitives/functions/Function1/makeFunction1s.C
    primitives/functions/Function1/reverseRamp/reverseRamp.C
    primitives/functions/Function1/linearRamp/linearRamp.C
    primitives/functions/Function1/quadraticRamp/quadraticRamp.C
    primitives/functions/Function1/quarterSineRamp/quarterSineRamp.C
    primitives/functions/Function1/quarterCosineRamp/quarterCosineRamp.C
    primitives/functions/Function1/halfCosineRamp/halfCosineRamp.C
    primitives/functions/Function1/exponentialSqrRamp/exponentialSqrRamp.C
    primitives/functions/Function1/normalise/normalise.C
    primitives/functions/Function1/squarePulse/squarePulse.C
    primitives/functions/Function1/Table/tableBase.C
    primitives/functions/Function1/Table/TableReader/makeTableReaders.C
    primitives/functions/Function1/unknownTypeFunction1/unknownTypeFunction1.C
    primitives/functions/Function1/omega1/omega1.C
    primitives/functions/Function1/Table/TableReader/Foam/tokenTupleN.C

    primitives/functions/Function2/Function2/Function2UnitConversions.C
    primitives/functions/Function2/makeFunction2s.C
)

set(strings primitives/strings)
target_sources(OpenFOAM PRIVATE
    ${strings}/string/string.C
    ${strings}/string/stringIO.C
    ${strings}/string/stringIOList.C
    ${strings}/verbatimString/verbatimString.C
    ${strings}/verbatimString/verbatimStringIO.C
    ${strings}/word/word.C
    ${strings}/word/wordIO.C
    ${strings}/word/wordIOList.C
    ${strings}/fileName/fileName.C
    ${strings}/fileName/fileNameIO.C
    ${strings}/variable/variable.C
    ${strings}/variable/variableIO.C
    ${strings}/functionName/functionName.C
    ${strings}/functionName/functionNameIO.C
    ${strings}/wordRe/wordRe.C
    ${strings}/lists/hashedWordList.C
    ${strings}/stringOps/stringOps.C
)
set(ops primitives/ops)
target_sources(OpenFOAM PRIVATE
    ${ops}/flipOp.C
)

set(sha1 primitives/hashes/SHA1)
target_sources(OpenFOAM PRIVATE
    primitives/hashes/Hasher/Hasher.C
    ${sha1}/SHA1.C
    ${sha1}/SHA1Digest.C
    primitives/randomGenerator/randomGenerator.C
)

set(ranges primitives/ranges)
target_sources(OpenFOAM PRIVATE
    ${ranges}/labelRange/labelRange.C
    ${ranges}/labelRange/labelRanges.C
    ${ranges}/scalarRange/scalarRange.C
    ${ranges}/scalarRange/scalarRanges.C
)

set(spatialVectorAlgebra primitives/spatialVectorAlgebra)
target_sources(OpenFOAM PRIVATE
    ${spatialVectorAlgebra}/SpatialVector/spatialVector/spatialVector.C
    ${spatialVectorAlgebra}/SpatialTensor/spatialTensor/spatialTensor.C
    ${spatialVectorAlgebra}/CompactSpatialTensor/compactSpatialTensor/compactSpatialTensor.C

)
target_sources(OpenFOAM PRIVATE
    primitives/polynomialEqns/cubicEqn/cubicEqn.C
    primitives/polynomialEqns/quadraticEqn/quadraticEqn.C

    primitives/Barycentric/barycentric/barycentric.C
    primitives/Barycentric2D/barycentric2D/barycentric2D.C

    containers/HashTables/HashTable/HashTableCore.C
    containers/Lists/SortableList/ParSortableListName.C
    containers/Lists/PackedList/PackedListCore.C
    containers/Lists/PackedList/PackedBoolList.C
    containers/Lists/ListOps/ListOps.C
    containers/LinkedLists/linkTypes/SLListBase/SLListBase.C
    containers/LinkedLists/linkTypes/DLListBase/DLListBase.C
)

set(Streams db/IOstreams)
target_sources(OpenFOAM PRIVATE
    ${Streams}/token/tokenIO.C
)

set(IOstreams ${Streams}/IOstreams)
target_sources(OpenFOAM PRIVATE
    ${IOstreams}/IOstream.C
    ${IOstreams}/Istream.C
    ${IOstreams}/Ostream.C
)

set(Sstreams ${Streams}/Sstreams)
target_sources(OpenFOAM PRIVATE
    ${Sstreams}/ISstream.C
    ${Sstreams}/OSstream.C
    ${Sstreams}/SstreamsPrint.C
    ${Sstreams}/readHexLabel.C
    ${Sstreams}/prefixOSstream.C
)

set(gzstream ${Streams}/gzstream)
target_sources(OpenFOAM PRIVATE
    ${gzstream}/gzstream.C
)


set(Fstreams ${Streams}/Fstreams)
target_sources(OpenFOAM PRIVATE
    ${Fstreams}/IFstream.C
    ${Fstreams}/OFstream.C
    ${Fstreams}/masterOFstream.C
)

set(Tstreams ${Streams}/Tstreams)
target_sources(OpenFOAM PRIVATE
    ${Tstreams}/ITstream.C
    ${Tstreams}/OTstream.C
)

set(StringStreams ${Streams}/StringStreams)
target_sources(OpenFOAM PRIVATE
    ${StringStreams}/StringStreamsPrint.C
)

set(Pstreams ${Streams}/Pstreams)
target_sources(OpenFOAM PRIVATE
    ${Pstreams}/UIPstream.C
    ${Pstreams}/IPstream.C
    # ${Pstreams}/UPstream.C in global.Cver
    ${Pstreams}/UPstreamCommsStruct.C
    ${Pstreams}/Pstream.C
    ${Pstreams}/UOPstream.C
    ${Pstreams}/OPstream.C
    ${Pstreams}/PstreamBuffers.C
)

set(dictionary db/dictionary)
target_sources(OpenFOAM PRIVATE
    ${dictionary}/dictionary.C
    ${dictionary}/dictionaryIO.C

    ${dictionary}/keyType/keyType.C
)

set(entry ${dictionary}/entry)
target_sources(OpenFOAM PRIVATE
    ${entry}/entry.C
    ${entry}/entryIO.C
)

set(primitiveEntry ${dictionary}/primitiveEntry)
target_sources(OpenFOAM PRIVATE
    ${primitiveEntry}/primitiveEntry.C
    ${primitiveEntry}/primitiveEntryIO.C
)

set(dictionaryEntry ${dictionary}/dictionaryEntry)
target_sources(OpenFOAM PRIVATE
    ${dictionaryEntry}/dictionaryEntry.C
    ${dictionaryEntry}/dictionaryEntryIO.C
)

set(dictionaryListEntry ${dictionary}/dictionaryListEntry)
target_sources(OpenFOAM PRIVATE
    ${dictionaryListEntry}/dictionaryListEntry.C
    ${dictionaryListEntry}/dictionaryListEntryIO.C

    ${dictionary}/mergeDictionaries/mergeDictionaries.C
)

set(functionEntries ${dictionary}/functionEntries)
target_sources(OpenFOAM PRIVATE
    ${functionEntries}/functionEntry/functionEntry.C
    ${functionEntries}/exitEntry/exitEntry.C
    ${functionEntries}/printEntry/printEntry.C
    ${functionEntries}/negEntry/negEntry.C
    ${functionEntries}/codeIncludeEntry/codeIncludeEntry.C
    ${functionEntries}/calcEntry/calcEntry.C
    ${functionEntries}/streamEntry/streamEntry.C
    ${functionEntries}/codeStream/codeStream.C
    ${functionEntries}/codeDict/codeDict.C
    ${functionEntries}/codeBlock/codeBlockEntry.C
    ${functionEntries}/codeBlockStream/codeBlockStreamEntry.C
    ${functionEntries}/codeBlockDict/codeBlockDictEntry.C
    ${functionEntries}/endCodeBlock/endCodeBlockEntry.C
    ${functionEntries}/includeEntry/includeEntry.C
    ${functionEntries}/includeEtcEntry/includeEtcEntry.C
    ${functionEntries}/includeFuncEntry/includeFuncEntry.C
    ${functionEntries}/includeIfPresentEntry/includeIfPresentEntry.C
    ${functionEntries}/inputModeEntry/inputModeEntry.C
    ${functionEntries}/removeEntry/removeEntry.C
    ${functionEntries}/ifeqEntry/ifeqEntry.C
    ${functionEntries}/ifEntry/ifEntry.C
    ${functionEntries}/elseEntry/elseEntry.C
    ${functionEntries}/elifEntry/elifEntry.C
    ${functionEntries}/endifEntry/endifEntry.C
)

set(IOdictionary db/IOobjects/IOdictionary)
target_sources(OpenFOAM PRIVATE
    ${IOdictionary}/IOdictionary.C
    ${IOdictionary}/IOdictionaryIO.C
    ${IOdictionary}/localIOdictionary.C
    ${IOdictionary}/timeIOdictionary.C
    ${IOdictionary}/systemDict.C
    db/IOobjects/decomposedBlockData/decomposedBlockData.C
)


set(IOobject db/IOobject)
target_sources(OpenFOAM PRIVATE
    ${IOobject}/IOobject.C
    ${IOobject}/IOobjectIO.C
    ${IOobject}/IOobjectReadHeader.C
    ${IOobject}/IOobjectWriteHeader.C
)

set(regIOobject db/regIOobject)
target_sources(OpenFOAM PRIVATE
    ${regIOobject}/regIOobject.C
    ${regIOobject}/regIOobjectRead.C
    ${regIOobject}/regIOobjectWrite.C
    db/IOobjectList/IOobjectList.C
    db/objectRegistry/objectRegistry.C
    db/CallbackRegistry/CallbackRegistryName.C
)


set(dll db/dynamicLibrary)
target_sources(OpenFOAM PRIVATE
    ${dll}/dlLibraryTable/dlLibraryTable.C
    ${dll}/dynamicCode/dynamicCode.C
    ${dll}/codedBase/codedBase.C
    ${dll}/compileTemplate/compileTemplate.C

    db/functionObjects/functionObject/functionObject.C
    db/functionObjects/functionObjectList/functionObjectList.C
    db/functionObjects/functionObjectList/wordAndDictionary.C
    db/functionObjects/writeFile/writeFile.C
    db/functionObjects/logFiles/logFiles.C
    db/functionObjects/logFile/logFile.C
    db/functionObjects/writeObjectsBase/writeObjectsBase.C
    db/functionObjects/writeLocalObjects/writeLocalObjects.C
    db/functionObjects/timeControl/timeControl.C
    db/functionObjects/timeControl/timeControlFunctionObject.C
    db/functionObjects/objectRegistryFunctionObject/objectRegistryFunctionObject.C
    db/functionObjects/regionFunctionObject/regionFunctionObject.C
)


set(Time db/Time)
target_sources(OpenFOAM PRIVATE
    ${Time}/controlIOdictionary/controlIOdictionary.C
    ${Time}/TimePaths.C
    ${Time}/TimeState.C
    ${Time}/Time.C
    ${Time}/TimeIO.C
    ${Time}/findTimes.C
    ${Time}/subCycleTime.C
    ${Time}/timeSelector.C

    ${Time}/instant/instant.C
)

set(userTime ${Time}/userTime)
target_sources(OpenFOAM PRIVATE
    ${userTime}/userTime/userTime.C
    ${userTime}/userTime/userTimeNew.C
    ${userTime}/realTime/realTime.C
    ${userTime}/engineTime/engineTime.C

    symbols/symbols.C
    dimensionSet/dimensionSet.C
    dimensionSet/dimensionSetIO.C
    dimensionSet/dimensionSets.C
    unitConversion/unitConversion.C
    unitConversion/unitConversionIO.C
    unitConversion/unitConversions.C
    unitConversion/namedUnitConversion.C
    unitConversion/namedUnitConversionIO.C
    dimensionedTypes/dimensionedScalar/dimensionedScalar.C
    dimensionedTypes/dimensionedSphericalTensor/dimensionedSphericalTensor.C
    dimensionedTypes/dimensionedSymmTensor/dimensionedSymmTensor.C
    dimensionedTypes/dimensionedTensor/dimensionedTensor.C

    matrices/solution/solution.C
)

set(scalarMatrices matrices/scalarMatrices)
target_sources(OpenFOAM PRIVATE
    ${scalarMatrices}/scalarMatrices.C
    ${scalarMatrices}/SVD/SVD.C
    ${scalarMatrices}/eigendecomposition/eigendecomposition.C
)

set(LUscalarMatrix matrices/LUscalarMatrix)
target_sources(OpenFOAM PRIVATE
    ${LUscalarMatrix}/LUscalarMatrix.C
    ${LUscalarMatrix}/procLduMatrix.C
    ${LUscalarMatrix}/procLduInterface.C
)


set(lduMatrix matrices/lduMatrix)
target_sources(OpenFOAM PRIVATE
    ${lduMatrix}/lduMatrix/lduMatrix.C
    ${lduMatrix}/lduMatrix/lduMatrixOperations.C
    ${lduMatrix}/lduMatrix/lduMatrixATmul.C
    ${lduMatrix}/lduMatrix/lduMatrixUpdateMatrixInterfaces.C
    ${lduMatrix}/lduMatrix/lduMatrixSolver.C
    ${lduMatrix}/lduMatrix/lduMatrixSmoother.C
    ${lduMatrix}/lduMatrix/lduMatrixPreconditioner.C

    ${lduMatrix}/solvers/diagonalSolver/diagonalSolver.C
    ${lduMatrix}/solvers/smoothSolver/smoothSolver.C
    ${lduMatrix}/solvers/PCG/PCG.C
    ${lduMatrix}/solvers/PBiCG/PBiCG.C
    ${lduMatrix}/solvers/PBiCGStab/PBiCGStab.C

    ${lduMatrix}/smoothers/noSmoother/noSmoother.C
    ${lduMatrix}/smoothers/GaussSeidel/GaussSeidelSmoother.C
    ${lduMatrix}/smoothers/symGaussSeidel/symGaussSeidelSmoother.C
    ${lduMatrix}/smoothers/nonBlockingGaussSeidel/nonBlockingGaussSeidelSmoother.C
    ${lduMatrix}/smoothers/DIC/DICSmoother.C
    ${lduMatrix}/smoothers/FDIC/FDICSmoother.C
    ${lduMatrix}/smoothers/DICGaussSeidel/DICGaussSeidelSmoother.C
    ${lduMatrix}/smoothers/DILU/DILUSmoother.C
    ${lduMatrix}/smoothers/DILUGaussSeidel/DILUGaussSeidelSmoother.C

    ${lduMatrix}/preconditioners/noPreconditioner/noPreconditioner.C
    ${lduMatrix}/preconditioners/diagonalPreconditioner/diagonalPreconditioner.C
    ${lduMatrix}/preconditioners/DICPreconditioner/DICPreconditioner.C
    ${lduMatrix}/preconditioners/FDICPreconditioner/FDICPreconditioner.C
    ${lduMatrix}/preconditioners/DILUPreconditioner/DILUPreconditioner.C
    ${lduMatrix}/preconditioners/GAMGPreconditioner/GAMGPreconditioner.C
)

set(lduAddressing ${lduMatrix}/lduAddressing)
target_sources(OpenFOAM PRIVATE
    ${lduAddressing}/lduAddressing.C
    ${lduAddressing}/lduInterface/lduInterface.C
    ${lduAddressing}/lduInterface/processorLduInterface.C
    ${lduAddressing}/lduInterface/cyclicLduInterface.C
)

set(lduInterfaceFields ${lduAddressing}/lduInterfaceFields)
target_sources(OpenFOAM PRIVATE
    ${lduInterfaceFields}/lduInterfaceField/lduInterfaceField.C
    ${lduInterfaceFields}/processorLduInterfaceField/processorLduInterfaceField.C
    ${lduInterfaceFields}/cyclicLduInterfaceField/cyclicLduInterfaceField.C
)

set(GAMG ${lduMatrix}/solvers/GAMG)
target_sources(OpenFOAM PRIVATE
    ${GAMG}/GAMGSolver.C
    ${GAMG}/GAMGSolverAgglomerateMatrix.C
    ${GAMG}/GAMGSolverInterpolate.C
    ${GAMG}/GAMGSolverScale.C
    ${GAMG}/GAMGSolverSolve.C
)

set(GAMGInterfaces ${GAMG}/interfaces)
target_sources(OpenFOAM PRIVATE
    ${GAMGInterfaces}/GAMGInterface/GAMGInterface.C
    ${GAMGInterfaces}/GAMGInterface/GAMGInterfaceNew.C
    ${GAMGInterfaces}/processorGAMGInterface/processorGAMGInterface.C
    ${GAMGInterfaces}/processorCyclicGAMGInterface/processorCyclicGAMGInterface.C
    ${GAMGInterfaces}/cyclicGAMGInterface/cyclicGAMGInterface.C
)

set(GAMGInterfaceFields ${GAMG}/interfaceFields)
target_sources(OpenFOAM PRIVATE
    ${GAMGInterfaceFields}/GAMGInterfaceField/GAMGInterfaceField.C
    ${GAMGInterfaceFields}/GAMGInterfaceField/GAMGInterfaceFieldNew.C
    ${GAMGInterfaceFields}/processorGAMGInterfaceField/processorGAMGInterfaceField.C
    ${GAMGInterfaceFields}/processorCyclicGAMGInterfaceField/processorCyclicGAMGInterfaceField.C
    ${GAMGInterfaceFields}/cyclicGAMGInterfaceField/cyclicGAMGInterfaceField.C

)
set(GAMGAgglomerations ${GAMG}/GAMGAgglomerations)
set(GAMGAgglomeration ${GAMGAgglomerations}/GAMGAgglomeration)
target_sources(OpenFOAM PRIVATE
    ${GAMGAgglomeration}/GAMGAgglomeration.C
    ${GAMGAgglomeration}/GAMGAgglomerateLduAddressing.C
)

set(pairGAMGAgglomeration ${GAMGAgglomerations}/pairGAMGAgglomeration)
target_sources(OpenFOAM PRIVATE
    ${pairGAMGAgglomeration}/pairGAMGAgglomeration.C
    ${pairGAMGAgglomeration}/pairGAMGAgglomerate.C
)

set(algebraicPairGAMGAgglomeration ${GAMGAgglomerations}/algebraicPairGAMGAgglomeration)
target_sources(OpenFOAM PRIVATE
    ${algebraicPairGAMGAgglomeration}/algebraicPairGAMGAgglomeration.C
)

set(dummyAgglomeration ${GAMGAgglomerations}/dummyAgglomeration)
target_sources(OpenFOAM PRIVATE
    ${dummyAgglomeration}/dummyAgglomeration.C
)


set(GAMGProcAgglomerations ${GAMG}/GAMGProcAgglomerations)
target_sources(OpenFOAM PRIVATE
    ${GAMGProcAgglomerations}/GAMGProcAgglomeration/GAMGProcAgglomeration.C
    ${GAMGProcAgglomerations}/all/allGAMGProcAgglomeration.C
    ${GAMGProcAgglomerations}/manual/manualGAMGProcAgglomeration.C
    ${GAMGProcAgglomerations}/sequential/sequentialGAMGProcAgglomeration.C
    ${GAMGProcAgglomerations}/none/noneGAMGProcAgglomeration.C
    ${GAMGProcAgglomerations}/pair/pairGAMGProcAgglomeration.C
)


target_sources(OpenFOAM PRIVATE
    meshes/lduMesh/lduMesh.C
    meshes/lduMesh/lduPrimitiveMesh.C
)

set(LduMatrix matrices/LduMatrix)
target_sources(OpenFOAM PRIVATE
    ${LduMatrix}/LduMatrix/lduMatrices.C
    ${LduMatrix}/LduMatrix/solverPerformance.C
    ${LduMatrix}/LduMatrix/LduInterfaceField/LduInterfaceFields.C
    ${LduMatrix}/Smoothers/lduSmoothers.C
    ${LduMatrix}/Preconditioners/lduPreconditioners.C
    ${LduMatrix}/Solvers/lduSolvers.C
)

set(primitiveShapes meshes/primitiveShapes)
target_sources(OpenFOAM PRIVATE
    ${primitiveShapes}/line/line.C
    ${primitiveShapes}/plane/plane.C
    ${primitiveShapes}/triangle/intersection.C
    ${primitiveShapes}/objectHit/pointIndexHitList.C
    ${primitiveShapes}/objectHit/pointIndexHitIOList.C
)

set(meshShapes meshes/meshShapes)
target_sources(OpenFOAM PRIVATE
    ${meshShapes}/edge/edge.C
    ${meshShapes}/edge/edgeIOList.C
)

set(face ${meshShapes}/face)
target_sources(OpenFOAM PRIVATE
    ${face}/face.C
    ${face}/faceIntersection.C
    ${face}/faceContactSphere.C
    ${face}/faceAreaInContact.C
    ${face}/faceIOList.C
)

set(cell ${meshShapes}/cell)
target_sources(OpenFOAM PRIVATE
    ${cell}/cell.C
    ${cell}/oppositeCellFace.C
    ${cell}/cellIOList.C
)

set(tetCell ${meshShapes}/tetCell)
target_sources(OpenFOAM PRIVATE
    ${tetCell}/tetCell.C
)

set(cellModeller ${meshShapes}/cellModeller)
target_sources(OpenFOAM PRIVATE
    ${cellModeller}/cellModeller.C
)

set(cellModel ${meshShapes}/cellModel)
target_sources(OpenFOAM PRIVATE
    ${cellModel}/cellModel.C
    ${cellModel}/cellModelIO.C
)

set(cellShape ${meshShapes}/cellShape)
target_sources(OpenFOAM PRIVATE
    ${cellShape}/cellShape.C
    ${cellShape}/cellShapeEqual.C
    ${cellShape}/cellShapeIO.C
    ${cellShape}/cellShapeIOList.C
)

target_sources(OpenFOAM PRIVATE
    meshes/Identifiers/patch/patchIdentifier.C
    meshes/Identifiers/patch/coupleGroupIdentifier.C

    meshes/meshObjects/meshObjects.C

    meshes/zones/pointZones/pointZone.C
    meshes/zones/pointZones/pointZoneList.C

    meshes/zones/cellZones/cellZone.C
    meshes/zones/cellZones/cellZoneList.C

    meshes/zones/faceZones/faceZone.C
    meshes/zones/faceZones/faceZoneList.C

    meshes/zoneGeneration/zoneSet/zoneSet.C
    meshes/zoneGeneration/zoneGenerator/zoneGenerator.C
    meshes/zoneGeneration/lookup/lookup.C

    meshes/zoneGeneration/zoneGeneratorList/zoneGeneratorList.C
    meshes/zoneGeneration/zonesGenerator/zonesGenerator.C

    meshes/zoneGeneration/generatedZoneSet/generatedZoneSet.C
)

set(polyMesh meshes/polyMesh)
set(polyPatches ${polyMesh}/polyPatches)
target_sources(OpenFOAM PRIVATE
    ${polyPatches}/polyPatch/polyPatch.C
    ${polyPatches}/polyPatch/polyPatchNew.C
)

set(basicPolyPatches ${polyPatches}/basic)
target_sources(OpenFOAM PRIVATE
    ${basicPolyPatches}/coupled/coupledPolyPatch.C
)

set(constraintPolyPatches ${polyPatches}/constraint)
target_sources(OpenFOAM PRIVATE
    ${constraintPolyPatches}/cyclic/cyclicTransform.C
    ${constraintPolyPatches}/cyclic/cyclicPolyPatch.C
    ${constraintPolyPatches}/cyclicSlip/cyclicSlipPolyPatch.C
    ${constraintPolyPatches}/empty/emptyPolyPatch.C
    ${constraintPolyPatches}/processorCyclic/processorCyclicPolyPatch.C
    ${constraintPolyPatches}/processor/processorPolyPatch.C
    ${constraintPolyPatches}/symmetryPlane/symmetryPlanePolyPatch.C
    ${constraintPolyPatches}/symmetry/symmetryPolyPatch.C
    ${constraintPolyPatches}/wedge/wedgePolyPatch.C
    ${constraintPolyPatches}/internal/internalPolyPatch.C
)

set(derivedPolyPatches ${polyPatches}/derived)
target_sources(OpenFOAM PRIVATE
    ${derivedPolyPatches}/wall/wallPolyPatch.C
)

set(polyBoundaryMesh ${polyMesh}/polyBoundaryMesh)
target_sources(OpenFOAM PRIVATE
    ${polyBoundaryMesh}/polyBoundaryMesh.C
    ${polyBoundaryMesh}/polyBoundaryMeshEntries.C
)

set(globalMeshData ${polyMesh}/globalMeshData)
target_sources(OpenFOAM PRIVATE
    ${globalMeshData}/globalMeshData.C
    ${globalMeshData}/globalPoints.C
    ${globalMeshData}/globalIndex.C
    ${globalMeshData}/commSchedule.C
    ${globalMeshData}/processorTopology.C

    ${polyMesh}/syncTools/syncTools.C
    ${polyMesh}/polyMeshTetDecomposition/polyMeshTetDecomposition.C
    ${polyMesh}/polyMeshTetDecomposition/tetIndices.C

    ${polyMesh}/polyMesh.C
    ${polyMesh}/polyMeshFromShapeMesh.C
    ${polyMesh}/polyMeshIO.C
    ${polyMesh}/polyMeshInitMesh.C
    ${polyMesh}/polyMeshClear.C
    ${polyMesh}/polyMeshUpdate.C

    ${polyMesh}/meshObjects/Residuals/residuals.C
    ${polyMesh}/meshObjects/cpuLoad/cpuLoad.C
)

set(primitiveMesh meshes/primitiveMesh)
target_sources(OpenFOAM PRIVATE
    ${primitiveMesh}/primitiveMesh.C
    ${primitiveMesh}/primitiveMeshCellCells.C
    ${primitiveMesh}/primitiveMeshCellCentresAndVols.C
    ${primitiveMesh}/primitiveMeshCellEdges.C
    ${primitiveMesh}/primitiveMeshCells.C
    ${primitiveMesh}/primitiveMeshClear.C
    ${primitiveMesh}/primitiveMeshEdgeCells.C
    ${primitiveMesh}/primitiveMeshEdgeFaces.C
    ${primitiveMesh}/primitiveMeshEdges.C
    ${primitiveMesh}/primitiveMeshFaceCentresAndAreas.C
    ${primitiveMesh}/primitiveMeshPointCells.C
    ${primitiveMesh}/primitiveMeshPointFaces.C
    ${primitiveMesh}/primitiveMeshPointPoints.C
    ${primitiveMesh}/primitiveMeshCellPoints.C
    ${primitiveMesh}/primitiveMeshCalcCellShapes.C
)

set(primitivePatch ${primitiveMesh}/primitivePatch)
target_sources(OpenFOAM PRIVATE
    ${primitivePatch}/patchZones.C
    ${primitivePatch}/walkPatch.C
)

set(cellMatcher ${meshShapes}/cellMatcher)
target_sources(OpenFOAM PRIVATE
    ${cellMatcher}/cellMatcher.C
    ${cellMatcher}/hexMatcher.C
    ${cellMatcher}/wedgeMatcher.C
    ${cellMatcher}/prismMatcher.C
    ${cellMatcher}/pyrMatcher.C
    ${cellMatcher}/tetWedgeMatcher.C
    ${cellMatcher}/tetMatcher.C
    ${cellMatcher}/degenerateMatcher.C
)

set(polyTopoChangeMap ${polyMesh}/polyTopoChangeMap)
target_sources(OpenFOAM PRIVATE
    ${polyTopoChangeMap}/polyTopoChangeMap.C
    ${polyTopoChangeMap}/faceMapper/faceMapper.C
    ${polyTopoChangeMap}/cellMapper/cellMapper.C
    ${polyTopoChangeMap}/mapAddedPolyMesh.C
)

set(polyMeshMap ${polyMesh}/polyMeshMap)
target_sources(OpenFOAM PRIVATE
    ${polyMeshMap}/polyMeshMap.C
)

set(polyDistributionMap ${polyMesh}/polyDistributionMap)
target_sources(OpenFOAM PRIVATE
    ${polyDistributionMap}/distributionMap.C
    ${polyDistributionMap}/distributionMapBase.C
    ${polyDistributionMap}/polyDistributionMap.C
    ${polyDistributionMap}/IOdistributionMap.C
)

set(PrimitivePatch ${primitiveMesh}/PrimitivePatch)
target_sources(OpenFOAM PRIVATE
    ${PrimitivePatch}/PrimitivePatchName.C
)


set(pointMesh meshes/pointMesh)
target_sources(OpenFOAM PRIVATE
    ${pointMesh}/pointMesh.C
)

set(pointMeshMapper ${pointMesh}/pointMeshMapper)
target_sources(OpenFOAM PRIVATE
    ${pointMeshMapper}/pointMapper.C
    ${pointMeshMapper}/pointPatchMapper.C
)

set(pointPatches ${pointMesh}/pointPatches)
target_sources(OpenFOAM PRIVATE
    ${pointPatches}/pointPatch/pointPatch.C
    ${pointPatches}/facePointPatch/facePointPatch.C
    ${pointPatches}/facePointPatch/facePointPatchNew.C
)

set(basicPointPatches ${pointPatches}/basic)
target_sources(OpenFOAM PRIVATE
    ${basicPointPatches}/coupled/coupledPointPatch.C
)

set(constraintPointPatches ${pointPatches}/constraint)
target_sources(OpenFOAM PRIVATE
    ${constraintPointPatches}/cyclic/cyclicPointPatch.C
    ${constraintPointPatches}/cyclicSlip/cyclicSlipPointPatch.C
    ${constraintPointPatches}/empty/emptyPointPatch.C
    ${constraintPointPatches}/processor/processorPointPatch.C
    ${constraintPointPatches}/processorCyclic/processorCyclicPointPatch.C
    ${constraintPointPatches}/symmetryPlane/symmetryPlanePointPatch.C
    ${constraintPointPatches}/symmetry/symmetryPointPatch.C
    ${constraintPointPatches}/wedge/wedgePointPatch.C
    ${constraintPointPatches}/internal/internalPointPatch.C
)

set(derivedPointPatches ${pointPatches}/derived)
target_sources(OpenFOAM PRIVATE
    ${derivedPointPatches}/coupled/coupledFacePointPatch.C
    ${derivedPointPatches}/wall/wallPointPatch.C
)

set(pointBoundaryMesh ${pointMesh}/pointBoundaryMesh)
target_sources(OpenFOAM PRIVATE
    ${pointBoundaryMesh}/pointBoundaryMesh.C
)

target_sources(OpenFOAM PRIVATE
    meshes/boundBox/boundBox.C

    meshes/treeBoundBox/treeBoundBox.C
)

set(meshTools meshes/meshTools)
target_sources(OpenFOAM PRIVATE
    ${meshTools}/matchPoints.C
)

target_sources(OpenFOAM PRIVATE
    fields/UniformDimensionedFields/uniformDimensionedFields.C

    fields/UniformGeometricFields/uniformGeometricFields.C
)

set(Fields fields/Fields)
target_sources(OpenFOAM PRIVATE
    ${Fields}/fieldMappers/fieldMapper/fieldMapper.C
    ${Fields}/fieldMappers/forwardFieldMapper/forwardFieldMapper.C
    ${Fields}/fieldMappers/forwardOrAssignFieldMapper/forwardOrAssignFieldMapper.C
    ${Fields}/fieldMappers/generalFieldMapper/generalFieldMapper.C
    ${Fields}/fieldMappers/identityFieldMapper/identityFieldMapper.C
    ${Fields}/fieldMappers/reverseFieldMapper/reverseFieldMapper.C
    ${Fields}/fieldMappers/reverseInterpolativeFieldMapper/reverseInterpolativeFieldMapper.C
    ${Fields}/fieldMappers/setSizeFieldMapper/setSizeFieldMapper.C
    ${Fields}/fieldMappers/setSizeAndZeroFieldMapper/setSizeAndZeroFieldMapper.C
    ${Fields}/labelField/labelField.C
    ${Fields}/scalarField/scalarField.C
    ${Fields}/vectorField/vectorField.C
    ${Fields}/vector2DField/vector2DField.C
    ${Fields}/sphericalTensorField/sphericalTensorField.C
    ${Fields}/diagTensorField/diagTensorField.C
    ${Fields}/symmTensorField/symmTensorField.C
    ${Fields}/tensorField/tensorField.C
    ${Fields}/quaternionField/quaternionField.C
    ${Fields}/triadField/triadField.C
    ${Fields}/complexFields/complexFields.C

    ${Fields}/labelField/labelIOField.C
    ${Fields}/labelField/labelFieldIOField.C
    ${Fields}/scalarField/scalarIOField.C
    ${Fields}/scalarField/scalarFieldIOField.C
    ${Fields}/vectorField/vectorIOField.C
    ${Fields}/vectorField/vectorFieldIOField.C
    ${Fields}/vector2DField/vector2DIOField.C
    ${Fields}/vector2DField/vector2DFieldIOField.C
    ${Fields}/sphericalTensorField/sphericalTensorIOField.C
    ${Fields}/sphericalTensorField/sphericalTensorFieldIOField.C
    ${Fields}/diagTensorField/diagTensorIOField.C
    ${Fields}/diagTensorField/diagTensorFieldIOField.C
    ${Fields}/symmTensorField/symmTensorIOField.C
    ${Fields}/symmTensorField/symmTensorFieldIOField.C
    ${Fields}/tensorField/tensorIOField.C
    ${Fields}/tensorField/tensorFieldIOField.C
    ${Fields}/quaternionField/quaternionIOField.C
    ${Fields}/triadField/triadIOField.C
    ${Fields}/transformField/transformField.C
)

set(pointPatchFields fields/pointPatchFields)
target_sources(OpenFOAM PRIVATE
    ${pointPatchFields}/pointPatchField/pointPatchFields.C
)

set(basicPointPatchFields ${pointPatchFields}/basic)
target_sources(OpenFOAM PRIVATE
    ${basicPointPatchFields}/calculated/calculatedPointPatchFields.C
    ${basicPointPatchFields}/coupled/coupledPointPatchFields.C
    ${basicPointPatchFields}/value/valuePointPatchFields.C
    ${basicPointPatchFields}/fixedValue/fixedValuePointPatchFields.C
    ${basicPointPatchFields}/zeroGradient/zeroGradientPointPatchFields.C
)

set(constraintPointPatchFields ${pointPatchFields}/constraint)
target_sources(OpenFOAM PRIVATE
    ${constraintPointPatchFields}/cyclic/cyclicPointPatchFields.C
    ${constraintPointPatchFields}/cyclicSlip/cyclicSlipPointPatchFields.C
    ${constraintPointPatchFields}/empty/emptyPointPatchFields.C
    ${constraintPointPatchFields}/processor/processorPointPatchFields.C
    ${constraintPointPatchFields}/processorCyclic/processorCyclicPointPatchFields.C
    ${constraintPointPatchFields}/symmetryPlane/symmetryPlanePointPatchFields.C
    ${constraintPointPatchFields}/symmetry/symmetryPointPatchFields.C
    ${constraintPointPatchFields}/wedge/wedgePointPatchFields.C
    ${constraintPointPatchFields}/internal/internalPointPatchFields.C
)

set(derivedPointPatchFields ${pointPatchFields}/derived)
target_sources(OpenFOAM PRIVATE
    ${derivedPointPatchFields}/slip/slipPointPatchFields.C
    ${derivedPointPatchFields}/fixedNormalSlip/fixedNormalSlipPointPatchFields.C
    ${derivedPointPatchFields}/uniformFixedValue/uniformFixedValuePointPatchFields.C
    ${derivedPointPatchFields}/codedFixedValue/codedFixedValuePointPatchFields.C
)

target_sources(OpenFOAM PRIVATE
    fields/GeometricFields/pointFields/pointFields.C

    meshes/bandCompression/bandCompression.C
    meshes/preservePatchTypes/preservePatchTypes.C
)

set(interpolationWeights interpolations/interpolationWeights)
target_sources(OpenFOAM PRIVATE
    ${interpolationWeights}/interpolationWeights/interpolationWeights.C
    ${interpolationWeights}/linearInterpolationWeights/linearInterpolationWeights.C
    ${interpolationWeights}/stepInterpolationWeights/stepInterpolationWeights.C
    ${interpolationWeights}/splineInterpolationWeights/splineInterpolationWeights.C
)

target_sources(OpenFOAM PRIVATE
    algorithms/pointInCell/pointInCell.C
    algorithms/indexedOctree/indexedOctreeName.C
    algorithms/indexedOctree/treeDataCell.C
    algorithms/indexedOctree/volumeType.C
    algorithms/dynamicIndexedOctree/dynamicIndexedOctreeName.C
    algorithms/dynamicIndexedOctree/dynamicTreeDataPoint.C
    algorithms/polygonTriangulate/polygonTriangulate.C
    algorithms/solutionControl/solutionControl.C
)

target_sources(OpenFOAM PRIVATE
    distributions/distribution/distribution.C
    distributions/distribution/distributionNew.C
    distributions/unintegrable/unintegrable.C
    distributions/unintegrable/unintegrableForNonZeroQ.C
    distributions/exponential/exponential.C
    distributions/fixedValue/fixedValue.C
    distributions/multiFixedValue/multiFixedValue.C
    distributions/multiNormal/multiNormal.C
    distributions/normal/normal.C
    distributions/RosinRammler/RosinRammler.C
    distributions/standardNormal/standardNormal.C
    distributions/tabulatedCumulative/tabulatedCumulative.C
    distributions/tabulatedDensity/tabulatedDensity.C
    distributions/uniform/uniform.C

    geometry/splines/BSpline/BSpline.C
    geometry/splines/CatmullRomSpline/CatmullRomSpline.C
)


target_link_libraries(OpenFOAM PRIVATE
    OSspecific
    Pstream_best
    ZLIB::ZLIB

    wmakeParse
)

target_compile_definitions(OpenFOAM PRIVATE
    LIB_NAME=libOpenFOAM.so
    WM_OPTIONS=${WM_OPTIONS}
)

target_include_directories(OpenFOAM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};OpenFOAM" PARENT_SCOPE)
