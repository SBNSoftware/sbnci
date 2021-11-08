SBN CI
======

Written by Andrew Scarff with copious amounts of help from Vito Di Benedetto and Dom Brailsford. Addtions and edits from Chris Hilgenberg,

This package is designed to work alongside the `lar_ci` continuous integration package to validate the current simulation chain using relevant metrics and plots. It contains the following main directories:

- `ups`: This contains all the generic setup to define this as a UPS product, including experiment code dependencies.
- `sbnci/Common/Modules`: This contains experiment agnostic, LArSoft analysis modules used to create TTrees of the metrics we want to compare.
- `sbnci/PlottingScripts`: This contains all experiment agnostic ROOT macros for reading in the TTrees and making ROOT plots from them, and then comparing those plots to reference plots.
- `sbnci/Common/scripts`:This contains the experiment agnostic bash scripts that are sourced by `lar_ci` to run the ROOT macros.
- `sbnci/<SBND/ICARUS>/Modules`: This contains experiment agnostic, LArSoft analysis modules used to create TTrees of the metrics we want to compare.
- `sbnci/<SBND/ICARUS>/PlottingScripts`: This contains all experiment specific ROOT macros for reading in the TTrees and making ROOT plots from them, and then comparing those plots to reference plots.
- `sbnci/<SBND/ICARUS>/JobConfigurations`: Experiment specific FHiCL files that configure modules in `sbnci/Common/Modules`
- `sbnci/<SBND/ICARUS>/thresholds`: This contains experiement specific text files with the threshold values for the chi2 comparisons made for each plot. The threshold values are used to determine the red-yellow-white colouring of the plots.

For more information see the Wiki page [here](https://sbnsoftware.github.io/sbn/sbnci_wiki/CI_Validation.html).
