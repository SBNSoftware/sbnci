#!/bin/bash

## Introduce environment variables specific to this validation
export ref_sbndcode_ana_hist=${ref_showervalidation_hist}
export normalise_plots=${normalise_plots_showervalidation}

## Name of plotting script. This should be the only line that needs changing for other CI chains.
plotScript="$SBNCI_DIR/scripts/sbnciplot_showervalidation.C"

## Passing the plotting script name and the input file name to sbnciplots.sh
source sbnciplots.sh $plotScript ${1}
