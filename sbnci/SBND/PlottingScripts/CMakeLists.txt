set( ${product}_script_dir "SBND/scripts" ) ### MIGRATE-CHANGED (migrate-3.20.00) - ${product}/+${version}/* -> ""
FILE( GLOB plot_scripts sbnciplot*.C  )
install( FILES ${plot_scripts} DESTINATION ${${product}_script_dir} )
