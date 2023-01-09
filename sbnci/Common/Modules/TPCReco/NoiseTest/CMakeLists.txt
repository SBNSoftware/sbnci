#simple_plugin(NoiseTest "module"
set(  MODULE_LIBRARIES
         larcore::Geometry_Geometry_service
         lardataobj::Simulation
         lardata::Utilities
         lardataobj::RawData
         lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
         sbndcode::Utilities_SignalShapingServiceSBND_service
         nusimdata::SimulationBase
         nurandom::RandomUtils_NuRandomService_service
         art_root_io::TFileService_service
         art::Framework_Services_Optional_RandomNumberGenerator_service
         art::Framework_Core
         art::Framework_Principal
         art::Framework_Services_Registry
         canvas::canvas
         messagefacility::MF_MessageLogger
         fhiclcpp::fhiclcpp
         cetlib::cetlib
         cetlib_except::cetlib_except
         CLHEP::CLHEP
         ROOT::Core
         art_root_io::tfile_support
         art::Framework_Services_Registry
         larcorealg::Geometry
         larsim::MCCheater_BackTrackerService_service
)

cet_build_plugin(NoiseTest art::Module SOURCE NoiseTest_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()

