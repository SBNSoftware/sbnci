#!/bin/bash

plotScript=`find ${FW_SEARCH_PATH//:/ } -maxdepth 1 -name sbnciplot_showervalidation.C -print -quit`

root -l -b -q $plotScript\(\"${1}\"\)

comparisonScript=`find ${FW_SEARCH_PATH//:/ } -maxdepth 1 -name CompareDataDistributions.C -print -quit`

root -l -b -q $comparisonScript\(\"${SBNDCODE_VERSION}\",\"${ref_sbndcode_version}\"\)
