
set(
  #EXCLUDE
  #  "LightenSimPhotons_module.cc"
  MODULE_LIBRARIES
    larreco::RecoAlg
    larsim::Simulation
    lardataobj::Simulation
    sbnci::Common_Modules_MCRecoUtils
    sbncode::OpDet_PDMapAlgSimple_tool
    lardataobj::RecoBase
    ROOT::Core
    ROOT::Physics
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    art_root_io::tfile_support ROOT::Core
    messagefacility::MF_MessageLogger
)

cet_build_plugin(LightenSimPhotons art::Module SOURCE LightenSimPhotons_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PDSRecoValidation art::Module SOURCE PDSRecoValidation_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PDSSimValidation  art::Module SOURCE PDSSimValidation_module.cc LIBRARIES ${MODULE_LIBRARIES})



install_fhicl()
install_source()
