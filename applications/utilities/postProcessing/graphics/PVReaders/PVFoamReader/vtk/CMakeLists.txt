VTK_MODULE_ADD_MODULE(PVFoamReader_VTK CLASSES vtkPVFoamReader)

#set_target_properties(
#    PVFoamReader_VTK
#    PROPERTIES
#    LIBRARY_OUTPUT_DIRECTORY "$ENV{PV_PLUGIN_PATH}"
#)

if(NOT PARAVIEW_QT_MAJOR_VERSION)
    set(PARAVIEW_QT_MAJOR_VERSION 5)
endif()

target_link_libraries(
    PVFoamReader_VTK
    PUBLIC
    OpenFOAM_Defines
    OSspecific
    OpenFOAM
    vtkPVFoam
    finiteVolume
    Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
    Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
    PRIVATE
    ParaView::pqApplicationComponents
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};PVFoamReader_VTK" PARENT_SCOPE)