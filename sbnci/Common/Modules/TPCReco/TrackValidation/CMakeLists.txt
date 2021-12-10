
art_make(MODULE_LIBRARIES
          larcorealg_Geometry
          larcore_Geometry_Geometry_service
          larsim_MCCheater_BackTrackerService_service
          larsim_MCCheater_ParticleInventoryService_service
          lardata_Utilities
          larreco_RecoAlg
          larreco_Calorimetry
          lardataalg_DetectorInfo
          nug4_ParticleNavigation
          sbndcode_RecoUtils
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ROOT::Core
          ROOT::Physics
          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
          ${MF_MESSAGELOGGER}
          ${ART_ROOT_IO_TFILE_SUPPORT}

        )

install_headers()
install_fhicl()
install_source()
