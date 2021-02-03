#!/bin/bash

## This script should be called from a specific .sh file for each validation process, e.g. showervalidation.sh

## Setting the plotting script location. 
plotScript=${1}

## Setting up reference file for comparison
echo "ref_sbndcode_ana_hist: ${ref_sbndcode_ana_hist}"
export ref_sbndcode_ana_hist_xrootd_path=$(echo ${ref_sbndcode_ana_hist} | sed 's#/pnfs/#root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/#g')
echo "ref_sbndcode_ana_hist_xrootd_path: ${ref_sbndcode_ana_hist_xrootd_path}"

echo "ref_sbndcode_version: ${ref_sbndcode_version}"
echo xrdcp -s ${ref_sbndcode_ana_hist_xrootd_path} $(basename ${ref_sbndcode_ana_hist_xrootd_path})
xrdcp -s ${ref_sbndcode_ana_hist_xrootd_path} $(basename ${ref_sbndcode_ana_hist_xrootd_path})
xrdcp_exitcode=$?
echo "xrdcp_exitcode: ${xrdcp_exitcode}"

mv -v ci_validation_histos.root ci_validation_histos_${ref_sbndcode_version}.root


## Run plotting script, producing ci_validation_histos.root
root -l -b -q $plotScript\(\"${2}\"\)


## Setup and run the comparison script.
comparisonScript="$SBNCI_DIR/scripts/CompareDataDistributions.C"

ln -v ci_validation_histos.root ci_validation_histos_${SBNDCODE_VERSION}.root

ls -lh

root -l -b -q $comparisonScript\(\"${SBNDCODE_VERSION}\",\"${ref_sbndcode_version}\"\)



## parse text files with validation results and create tables out of them
#  for now we only have ComparisonChiSquare.txt
summary_table_list="ComparisonChiSquare.txt"

for summary_table in ${summary_table_list}; do
    summary_table_noext=${summary_table//.*} # summary_table w/o file extension

    ## table format is as follow
    # first row is the table header
    # other rows are for data
    # element separator is semicolon
    # default cell color background is white
    # to have custom cell background we can append HTML colors to the cell value separated by #
    # a typical element in the row is like
    # val 1; val 2; val 3#COLOR

    echo "Sbndcode Release;Parameter;Chi2" | tee ./${summary_table_noext}_html_table
    ## for Chi2 values we define the background as follow:
    #  <0.25 = green, <0.5 = orange and >=0.5 = red
    awk 'BEGIN{GREEN="#008000"; ORANGE="#e79903"; RED="#FF0000"}
        {if ($2<0.25) COLOR=GREEN ; else if ($2<0.5) COLOR=ORANGE; else COLOR=RED;
        { print "'"${ref_sbndcode_version} vs ${SBNDCODE_VERSION}"';"$1";"$2COLOR}}' ./${summary_table} |
    tee -a ./${summary_table_noext}_html_table

done
