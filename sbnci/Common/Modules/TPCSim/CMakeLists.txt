
set(MODULE_LIBRARIES
          larreco::RecoAlg
          lardataobj::RecoBase
          art::Framework_Services_Registry
          ROOT::Core
          ROOT::Physics
          art_root_io::TFileService_service
          messagefacility::MF_MessageLogger
          art_root_io::tfile_support
          nusimdata::SimulationBase
          sbnci::Common_Modules_MCRecoUtils
        )

cet_build_plugin(TPCSimValidation art::Module SOURCE TPCSimValidation_module.cc LIBRARIES ${MODULE_LIBRARIES} )

install_fhicl()
install_source()
