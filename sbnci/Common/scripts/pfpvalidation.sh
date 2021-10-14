#!/bin/bash

# Set code name and version to have experiment agnostic env vars
source sbnci_setcodename.sh

## Introduce environment variables specific to this validation
## These are defined in the corresponding lar_ci cfg file
export ref_ana_hist=${ref_pfp_hist}
export normalise_plots=${normalise_plots_pfp}

## Name of plotting script. This should be the only line that needs changing for other CI chains.
plotScript="${expCIDir}/scripts/sbnciplot_pfpvalidation.C"

## Passing the plotting script name and the input file name to sbnciplots.sh
source sbnciplots.sh $plotScript ${1}
