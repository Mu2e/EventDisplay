# the code branches on the version, so have to provide in the right format
string(REPLACE "." ";" version_list "${ROOT_VERSION}")
list(GET version_list 0 a)
list(GET version_list 1 b)
list(GET version_list 2 c)
math(EXPR NRV "${a}*1000000 + ${b}*10000 + ${c}*100")
message("---> Converted ROOT_VERSION ${ROOT_VERSION} to NUMERIC_ROOT_VERSION ${NRV} ")
add_compile_options("-DNUMERIC_ROOT_VERSION=${NRV}")

cet_make_library(
    SOURCE
      CollectionFiller.cc
      DataInterface.cc
      TrackerCalo2DViews.cc
      EventDisplayManager.cc
      GUI.cc
      MCInterface.cc
      MainWindow.cc
      PrintInfo.cc
      TextSelect.cc
    LIBRARIES PUBLIC
      KinKal_Trajectory
      Offline::CRVReco
      Offline::CRVResponse
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::DataProducts
      Offline::KinKalGeom
      Offline::Mu2eKinKal
      Offline::TrackerConditions
      Offline::Mu2eUtilities
      Offline::GeneralUtilities
      ROOT::Eve
      ROOT::ROOTEve
      ROOT::ROOTWebDisplay

    LIBRARY_NAME eventdisplay
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/classes.h
     DICTIONARY_LIBRARIES
      event-display::eventdisplay
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::DataProducts
      Offline::Mu2eUtilities
      Offline::GeneralUtilities
)

cet_build_plugin(Mu2eEventDisplay art::module
    REG_SOURCE   Mu2eEventDisplay_module.cc
    LIBRARIES REG
      event-display::eventdisplay
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::DataProducts
      Offline::Mu2eUtilities
      Offline::GeneralUtilities
      Offline::RecoDataProducts
      art_root_io::TFileService_service
      ROOT::Postscript
)

install_source(SUBDIRS src)
