#!/bin/bash

##THISCIDIR=$(eval "echo \${${PROJ_PREFIX}_CI_DIR}")

plotScript=`find ${FW_SEARCH_PATH//:/ } -maxdepth 1 -name sbnciplot_showervalidation.C -print -quit`

root -l -b -q $plotScript

##root -l -b -q ${THISCIDIR}/test/CompareDataDistributions.C\(\"${DUNETPC_VERSION}\",\"${ref_dunetpc_version}\"\)
