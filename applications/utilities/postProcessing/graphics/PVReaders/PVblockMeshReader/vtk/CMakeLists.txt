VTK_MODULE_ADD_MODULE(PVblockMeshReader_VTK CLASSES vtkPVblockMeshReader)

#set_target_properties(
#    PVblockMeshReader_VTK
#    PROPERTIES
#    LIBRARY_OUTPUT_DIRECTORY "$ENV{PV_PLUGIN_PATH}"
#)

if(NOT PARAVIEW_QT_MAJOR_VERSION)
    set(PARAVIEW_QT_MAJOR_VERSION 5)
endif()

target_link_libraries(PVblockMeshReader_VTK PUBLIC
    OpenFOAM
    OpenFOAM_Defines
    OSspecific

    vtkPVblockMesh
    ParaView::pqApplicationComponents
    Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
    Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
    LINK_PRIVATE
    ParaView::pqApplicationComponents
)
