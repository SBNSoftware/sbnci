
set(MODULE_LIBRARIES
          larcorealg::Geometry
          larcore::Geometry_Geometry_service
          larsim::MCCheater_BackTrackerService_service
          larsim::MCCheater_ParticleInventoryService_service
          lardata::Utilities
          larreco::RecoAlg
          larreco::Calorimetry
          lardataalg::DetectorInfo
          nug4::ParticleNavigation
          sbndcode_RecoUtils
          ROOT::Core
          ROOT::Physics
          nusimdata::SimulationBase
          messagefacility::MF_MessageLogger
          art::Framework_Services_Registry
          art_root_io::tfile_support
          art_root_io::TFileService_service
        )

cet_build_plugin(TrackValidation art::Module SOURCE TrackValidation_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
