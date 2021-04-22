Modules
=======

This directory contains the LArSoft analysis modules used to create metrics in `TTrees` for physics validation. 
It also contains utilities to be used in the analysis modules.

Each analysis module should be accompanied by a prolog-style fcl file defining a table with the configuration
for that module. The fcls for running these modules during a validation job should live in `JobConfigurations`.
