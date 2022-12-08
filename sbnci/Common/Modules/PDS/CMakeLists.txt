
set(
  #EXCLUDE
  #  "LightenSimPhotons_module.cc"
  MODULE_LIBRARIES
    larreco::RecoAlg
    larsim::Simulation
    lardataobj::Simulation
    sbnci_Common_Modules_MCRecoUtils
    lardataobj::RecoBase
    ROOT::Core
    ROOT::Physics
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    art_root_io::tfile_support ROOT::Core
    messagefacility::MF_MessageLogger
)

cet_build_plugin(LightenSimPhotons art::Module SOURCE LightenSimPhotons_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
