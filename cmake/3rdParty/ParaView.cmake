if(${ParaView_TYPE} STREQUAL "none")
    message(STATUS "Skip finding openfoam because ParaView_TYPE = ${ParaView_TYPE}")
    return()
endif ()
#find_package(VTK CONFIG REQUIRED)
find_package(ParaView 5.7 CONFIG REQUIRED)

# ParaView will import vtk; while vtk imports tons of deps. adios2sys_objects is handled with freaking bugs.
# Add dummy target to make vtk happy.
if (NOT TARGET adios2sys_objects)
    add_library(adios2sys_objects OBJECT
        ${CMAKE_CURRENT_LIST_DIR}/dummy.cpp
    )
endif ()