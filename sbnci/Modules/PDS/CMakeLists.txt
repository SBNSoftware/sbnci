
art_make(
  EXCLUDE
    "LightenSimPhotons_module.cc"
  MODULE_LIBRARIES
    larreco_RecoAlg
    sbnci_Modules_MCRecoUtils
          lardataobj_RecoBase
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ROOT::Core
          ROOT::Physics
          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
          ${MF_MESSAGELOGGER}
    ${ART_ROOT_IO_TFILE_SUPPORT}
        )

simple_plugin(LightenSimPhotons module
  larreco_RecoAlg
  sbnci_Modules_MCRecoUtils
        lardataobj_RecoBase
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
