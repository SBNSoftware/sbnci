#set( ${product}_script_dir "${product}/${version}/Common/scripts" )
set( ${product}_script_dir "${sbnci_UPS_PRODUCT_NAME}/${sbnci_UPS_PRODUCT_VERSION}/Common/scripts" )
FILE( GLOB plot_scripts CompareDataDistributions.C )
install( FILES ${plot_scripts} DESTINATION ${${product}_script_dir} )
