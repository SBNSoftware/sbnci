#!/bin/bash

## This script should be called from a specific .sh file for each validation process, e.g. showervalidation.sh

## Setting the plotting script location. 
plotScript=${1}

## Setting up reference file for comparison
echo "ref_ana_hist: ${ref_ana_hist}"
export ref_ana_hist_xrootd_path=$(echo ${ref_ana_hist} | sed 's#/pnfs/#root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/#g')
echo "ref_ana_hist_xrootd_path: ${ref_ana_hist_xrootd_path}"
echo "normalise_plots: ${normalise_plots}"

ref_version=""
if [ ${codeName} == "sbndcode" ]
then
  ref_version="${ref_sbndcode_version}"
elif [ ${codeName} == "icaruscode" ]
then
  ref_version="${ref_icaruscode_version}"
fi

echo "ref_${codeName}_version: ${ref_version}"
echo xrdcp -s ${ref_ana_hist_xrootd_path} $(basename ${ref_ana_hist_xrootd_path})
xrdcp -s ${ref_ana_hist_xrootd_path} $(basename ${ref_ana_hist_xrootd_path})
xrdcp_exitcode=$?
echo "xrdcp_exitcode: ${xrdcp_exitcode}"

mv -v ci_validation_histos.root ci_validation_histos_ref_${ref_version}.root


## Run plotting script, producing ci_validation_histos.root
root -l -b -q $plotScript\(\"${2}\"\)


## Setup and run the comparison script.
#comparisonScript="$SBNCI_DIR/Common/PlottingScripts/CompareDataDistributions.C"
comparisonScript="$SBNCI_DIR/scripts/CompareDataDistributions.C"

ln -v ci_validation_histos.root ci_validation_histos_${codeVersion}.root

ls -lh

root -l -b -q $comparisonScript\(\"${codeVersion}\",\"${ref_version}\",${normalise_plots:-false}\)
