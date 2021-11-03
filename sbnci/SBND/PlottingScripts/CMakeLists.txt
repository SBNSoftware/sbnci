set( ${product}_script_dir "${product}/${version}/scripts" )
FILE( GLOB plot_scripts sbnciplot*.C  )
install( FILES ${plot_scripts} DESTINATION ${${product}_script_dir} )
