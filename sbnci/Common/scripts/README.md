scripts
=======

This directory contains the bash scripts used to run validation jobs. The job specific scripts
set the environmental variables specific to their validation process and then call `sbnciplots.sh`.

`sbnciplots.sh` is the high level script that calls the plotting and comparison scripts.