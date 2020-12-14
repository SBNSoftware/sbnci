#!/bin/bash

plotScript="$SBNCI_DIR/scripts/sbnciplot_showervalidation.C"

root -l -b -q $plotScript\(\"${1}\"\)

comparisonScript="$SBNCI_DIR/scripts/CompareDataDistributions.C"

ln -v ci_validation_histos.root ci_validation_histos_${SBNDCODE_VERSION}.root

ls -lh

root -l -b -q $comparisonScript\(\"${SBNDCODE_VERSION}\",\"${ref_sbndcode_version}\"\)
