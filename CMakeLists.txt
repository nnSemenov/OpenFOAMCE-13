cmake_minimum_required(VERSION 3.28)
project(Mikeno
    VERSION 14
    LANGUAGES CXX
    DESCRIPTION "Mikeno is a fork of OpenFOAM, Frankensteined for chemical engineering usage."
    HOMEPAGE_URL "https://github.com/nnSemenov/Mikeno"
)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(WM_PRECISION_OPTION DP CACHE STRING "Precision of float point number. SP | DP")

set(WM_LABEL_SIZE 32 CACHE STRING "Size of label. 32 | 64")

set(WM_MPLIB SYSTEMOPENMPI CACHE STRING "MPI implementation. SYSTEMOPENMPI | OPENMPI | SYSTEMMPI | MPICH
    | MPICH-GM | HPMPI | MPI | FJMPI | QSMPI | SGIMPI | INTELMPI")

set(WM_ARCH linux64 CACHE STRING "Computer arch")
set(WM_ARCH_OPTION 64 CACHE STRING "Computer arch option")
set(WM_OSTYPE posix CACHE STRING "OS type")

#- Decomposition libraries:
#    <library>_TYPE = none | system | ThirdParty (OpenFOAM)
set(SCOTCH_TYPE system CACHE STRING "Scotch type")
set(METIS_TYPE system CACHE STRING "")
set(PARMETIS_TYPE system CACHE STRING "")
set(ZOLTAN_TYPE system CACHE STRING "")

#    ParaView_TYPE = none | system | paraviewopenfoam | ThirdParty (OpenFOAM)
set(ParaView_TYPE system CACHE STRING "")

#- Floating-point signal handling:
#    set or unset
set(FOAM_SIGFPE "" CACHE STRING "")
#export FOAM_SIGFPE=

#- Memory initialisation:
#    set or unset
set(FOAM_SETNAN "" CACHE STRING "")
#export FOAM_SETNAN=

option(MIKENO_UNIT_TEST "Enable unit test for Mikeno" ON)

include(cmake/3rdParty.cmake)

include(cmake/process_options.cmake)

include(cmake/process_Cver.cmake)
include(cmake/lnInclude.cmake)

# libraries that are installed to lib.
set(FOAM_regular_libraries)
# libraries that don't install to lib, but different destination
set(FOAM_special_libraries)
# executables installed to bin
set(FOAM_regular_executables)

add_subdirectory(src)
add_subdirectory(applications)

include(cmake/install.cmake)