set(MODULE_LIBRARIES
        larsim::Utils
        nusimdata::SimulationBase
        lardataobj::RecoBase
        sbnci::Common_Modules_MCRecoUtils
        nusimdata::SimulationBase
        art::Framework_Services_Registry
        art_root_io::tfile_support
        art_root_io::TFileService_service
        ROOT::Core
        ROOT::Tree
        sbndcode::GeoWrappers
        )

cet_build_plugin(RecoEff art::Module SOURCE RecoEff_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_fhicl()
install_source()
