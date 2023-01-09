art_make_library(
             SOURCE
                 RecoUtils.cc
                 ShowerUtils.cc
                 TrackUtils.cc
             LIBRARY_NAME
                 sbnci::Common_Modules_MCRecoUtils
             LIBRARIES
                 larcorealg::Geometry
                 larcore::Geometry_Geometry_service
                 larsim::Simulation lardataobj::Simulation
                 larsim::MCCheater_BackTrackerService_service
                 larsim::MCCheater_ParticleInventoryService_service
                 lardata::Utilities
                 lardataobj::RawData
                 lardataobj::RecoBase
                 lardata::RecoObjects
                 nusimdata::SimulationBase
                 cetlib::cetlib
                 cetlib_except::cetlib_except
                 #${ART_FRAMEWORK_SERVICES_OPTIONAL}
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support
                 art_root_io::TFileService_service
                 art::Persistency_Common
                 art::Persistency_Provenance
                 art::Utilities
                 canvas::canvas
                 messagefacility::MF_MessageLogger
                 fhiclcpp::fhiclcpp
                 ROOT::Geom
                 ROOT::XMLIO
                 ROOT::Gdml
                 ROOT::Core
                 #${Boost_SYSTEM_LIBRARY}


        )

install_headers()
install_fhicl()
install_source()
