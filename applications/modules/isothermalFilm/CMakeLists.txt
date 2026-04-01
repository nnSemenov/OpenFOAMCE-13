lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(filmCompressibleMomentumTransportModels)

set(target_name isothermalFilmSolver)

add_library(${target_name} SHARED
    patches/filmFvPatch/filmFvPatch.C

    patches/filmWall/filmWallPolyPatch/filmWallPolyPatch.C
    patches/filmWall/filmWallPointPatch/filmWallPointPatch.C
    patches/filmWall/filmWallFvPatch/filmWallFvPatch.C

    patches/filmSurface/filmSurfacePolyPatch/filmSurfacePolyPatch.C
    patches/filmSurface/filmSurfacePointPatch/filmSurfacePointPatch.C
    patches/filmSurface/filmSurfaceFvPatch/filmSurfaceFvPatch.C

    patches/mappedFilmWall/mappedFilmWallFvPatch/mappedFilmWallFvPatch.C
    patches/mappedFilmWall/mappedFilmWallPointPatch/mappedFilmWallPointPatch.C
    patches/mappedFilmWall/mappedFilmWallPolyPatch/mappedFilmWallPolyPatch.C

    patches/mappedFilmSurface/mappedFilmSurfacePolyPatch/mappedFilmSurfacePolyPatch.C
    patches/mappedFilmSurface/mappedFilmSurfacePointPatch/mappedFilmSurfacePointPatch.C
    patches/mappedFilmSurface/mappedFilmSurfaceFvPatch/mappedFilmSurfaceFvPatch.C

    derivedFvPatchFields/alphaOne/alphaOneFvPatchScalarField.C
    derivedFvPatchFields/filmContactAngle/filmContactAngleFvPatchScalarField.C
    derivedFvPatchFields/mappedFilmPressure/mappedFilmPressureFvPatchScalarField.C
    derivedFvPatchFields/filmSurfaceVelocity/filmSurfaceVelocityFvPatchVectorField.C

    filmGaussGrad/filmGaussGrads.C

    contactForce.C
    continuityPredictor.C
    momentumPredictor.C
    correctAlpha.C
    isothermalFilm.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    # ${LIB_SRC}/MomentumTransportModels/compressible/lnInclude
#    ${LIB_SRC}/twoPhaseModels/interfaceProperties/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    fluidThermophysicalModels
    momentumTransportModels
    compressibleMomentumTransportModels
    filmCompressibleMomentumTransportModels
    interfaceProperties
    finiteVolume
    meshTools
    sampling
    fvModels
    fvConstraints
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

add_subdirectory(fvModels/filmVoFTransfer)
add_subdirectory(fvModels/filmCloudTransfer)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


