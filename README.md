SBN CI
======

Written by Andrew Scarff with copious amounts of help from Vito Di Benedetto and Dom Brailsford.

This package is designed to work alongside the `lar_ci` continuous integration package to validate the current simulation chain using relevant metrics and plots. It contains the following main directories:

- `ups`: This contains all the generic setup to define this as a UPS product.
- `sbnci/Modules`: This contains LArSoft analysis modules used to create TTrees of the metrics we want to compare.
- `sbnci/PlottingScripts`: This contains all the ROOT macros for reading in the TTrees and making ROOT plots from them, and then comparing those plots to reference plots.
- `sbnci/scripts`: This contains the bash scripts that are sourced by `lar_ci` to run the ROOT macros.
- `sbnci.thresholds`: This contains text files with the threshold values for the chi2 comparisons made for each plot. The threshold values are used to determine the red-yellow-white colouring of the plots.

For more information see the Wiki page [here](https://sbnsoftware.github.io/sbndcode_wiki/CI_Validation.html).
