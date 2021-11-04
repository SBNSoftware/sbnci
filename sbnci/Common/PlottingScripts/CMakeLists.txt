set( ${product}_script_dir "${product}/${version}/Common/scripts" )
FILE( GLOB plot_scripts CompareDataDistributions.C )
install( FILES ${plot_scripts} DESTINATION ${${product}_script_dir} )
