
set(MODULE_LIBRARIES
    larreco::RecoAlg
    lardataobj::RecoBase
    #${ART_FRAMEWORK_SERVICES_REGISTRY}
    art::Framework_Services_Registry
    ROOT::Core
    ROOT::Physics
    #${ART_ROOT_IO_TFILESERVICE_SERVICE}
    art_root_io::TFileService_service
    #${MF_MESSAGELOGGER}
    messagefacility::MF_MessageLogger
    #${ART_ROOT_IO_TFILE_SUPPORT}
    art_root_io::tfile_support
    sbnobj::Common_CRT
    nusimdata::SimulationBase
    sbnci::Common_Modules_MCRecoUtils
)

cet_build_plugin(CRTValidation art::Module SOURCE CRTValidation_module.cc LIBRARIES ${MODULE_LIBRARIES})


#install_headers()
install_fhicl()
install_source()
