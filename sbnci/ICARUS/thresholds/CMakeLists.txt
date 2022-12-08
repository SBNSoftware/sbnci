set( ${product}_threshold_dir "ICARUS/thresholds" ) ### MIGRATE-CHANGED (migrate-3.20.00) - ${product}/+${version}/* -> ""
FILE( GLOB threshold_files *ComparisonChiSquareTh.txt )
install( FILES ${threshold_files} DESTINATION ${${product}_threshold_dir} )
