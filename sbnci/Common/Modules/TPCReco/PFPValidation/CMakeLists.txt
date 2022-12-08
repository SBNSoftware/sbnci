set(MODULE_LIBRARIES
          sbnci_Common_Modules_MCRecoUtils
          lardataobj::RecoBase
          larsim::Utils
	  lardataalg::DetectorInfo
          nusimdata::SimulationBase
          ROOT::Core
          ROOT::Physics
          ROOT::Tree
          messagefacility::MF_MessageLogger
          art::Framework_Services_Registry
          art_root_io::tfile_support
          art_root_io::TFileService_service
        )

cet_build_plugin(PFPSliceValidationCI art::Module SOURCE PFPSliceValidationCI_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PFPValidationCI art::Module SOURCE PFPValidationCI_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
