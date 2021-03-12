#!/bin/bash

## Allows for different naming structure in cfg file
export ref_sbndcode_ana_hist=$(echo ${ref_pfpslice_hist})

## Name of plotting script. This should be the only line that needs changing for other CI chains.
plotScript="$SBNCI_DIR/scripts/sbnciplot_pfpslicevalidation.C"

## Passing the plotting script name and the input file name to sbnciplots.sh
source sbnciplots.sh $plotScript ${1}
