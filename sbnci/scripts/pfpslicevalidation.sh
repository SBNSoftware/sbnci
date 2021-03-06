#!/bin/bash

## Introduce environment variables specific to this validation
export ref_sbndcode_ana_hist=${ref_pfpslice_hist}
export normalise_plots=${normalise_plots_pfpslice}

## Name of plotting script. This should be the only line that needs changing for other CI chains.
plotScript="$SBNCI_DIR/scripts/sbnciplot_pfpslicevalidation.C"

## Passing the plotting script name and the input file name to sbnciplots.sh
source sbnciplots.sh $plotScript ${1}
