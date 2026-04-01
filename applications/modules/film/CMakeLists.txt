lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(filmThermophysicalTransportModels)

set(target_name filmSolver)

add_library(${target_name} SHARED
    film.C
    thermophysicalPredictor.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE
    ${LIB_SRC}/physicalProperties/lnInclude
    ${LIB_SRC}/thermophysicalModels/basic/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    ${LIB_SRC}/ThermophysicalTransportModels/fluid/lnInclude

    ${FOAM_MODULES}/isothermalFilm/filmCompressibleMomentumTransportModels/lnInclude
)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific

    isothermalFilmSolver

    fluidThermophysicalTransportModel
    thermophysicalTransportModel
    filmThermophysicalTransportModels
    coupledThermophysicalTransportModels
    momentumTransportModels
    compressibleMomentumTransportModels
    finiteVolume
    meshTools
    sampling
    fvModels
    fvConstraints
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


