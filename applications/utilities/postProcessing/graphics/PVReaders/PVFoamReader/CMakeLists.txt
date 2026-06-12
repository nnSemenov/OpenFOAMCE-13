project(PVFoamReader)

find_package(ParaView REQUIRED)

if(${ParaView_VERSION} VERSION_LESS 5.7)
    # ParaView-5.7.x and lower will prepend "lib" to the target name
    set(TARGET_NAME PVFoamReader_SM)
else()
    set(TARGET_NAME libPVFoamReader_SM)
endif()

if(${ParaView_VERSION} VERSION_LESS 5.7)
    # ParaView-5.7.x and lower requires the paraview include file
    include(${PARAVIEW_USE_FILE})
endif()


# Add the plugin
if(${ParaView_VERSION} VERSION_LESS 5.7)
    add_paraview_plugin(
        ${TARGET_NAME}
        "1.0"
        SERVER_MANAGER_XML PVFoamReader_SM.xml
        SERVER_MANAGER_SOURCES vtk/vtkPVFoamReader.cxx
    )
else()
    # Paraview-5.7.x and higher builds the vtk module separately
    paraview_add_plugin(
        ${TARGET_NAME}
        VERSION "1.0"
        SERVER_MANAGER_XML PVFoamReader_SM.xml
        MODULES PVFoamReader_VTK
        MODULE_FILES "vtk/vtk.module"
    )
endif()

# Set the output library destination to the plugin directory
#set_target_properties(
#    ${TARGET_NAME}
#    PROPERTIES
#    LIBRARY_OUTPUT_DIRECTORY "$ENV{PV_PLUGIN_PATH}"
#)


target_link_libraries(${TARGET_NAME} PUBLIC
    OpenFOAM_Defines
    OSspecific
    OpenFOAM

    vtkPVFoam
    finiteVolume
)

# Build the plugin
if(${ParaView_VERSION} VERSION_LESS 5.5)
    target_link_libraries(${TARGET_NAME}
        PRIVATE
        pqApplicationComponents
    )
elseif(${ParaView_VERSION} VERSION_LESS 5.7)
    # Assume QT version 5 if it hasn't been defined
    if(NOT PARAVIEW_QT_MAJOR_VERSION)
        set(PARAVIEW_QT_MAJOR_VERSION 5)
    endif()
    # Paraview-5.5.x and higher needs QT libraries listed in the link command
    # in order for the headers to be available
    target_link_libraries(${TARGET_NAME}
        PUBLIC
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
        PRIVATE
        pqApplicationComponents
    )
else()
    # Assume QT version 5 if it hasn't been defined
    if(NOT PARAVIEW_QT_MAJOR_VERSION)
        set(PARAVIEW_QT_MAJOR_VERSION 5)
    endif()
    # Paraview-5.5.x and higher needs QT libraries listed in the link command
    # in order for the headers to be available
    target_link_libraries(${TARGET_NAME}
        PUBLIC
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
        PRIVATE
        ParaView::pqApplicationComponents
    )
endif()

set(FOAM_regular_libraries "${FOAM_regular_libraries};${TARGET_NAME}" PARENT_SCOPE)
