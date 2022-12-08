
set(MODULE_LIBRARIES
    larreco::RecoAlg
    sbnci_Common_Modules_MCRecoUtils
    lardataobj::RecoBase
    nusimdata::SimulationBase
    ROOT::Core
    ROOT::Physics
    messagefacility::MF_MessageLogger
    art::Framework_Services_Registry
    art_root_io::tfile_support
    art_root_io::TFileService_service
)

cet_build_plugin(ShowerValidation art::Module SOURCE ShowerValidation_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
