
set( ${product}_threshold_dir "${product}/${version}/thresholds" )
FILE( GLOB threshold_files *ComparisonChiSquareTh.txt )
install( FILES ${threshold_files} DESTINATION ${${product}_threshold_dir} )
