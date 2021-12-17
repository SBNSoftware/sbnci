#!/bin/bash

source sbnci_setcodename.sh # figure out if we're doing SBND or ICARUS stuff
refdir=/pnfs/${expName}/persistent/ContinuousIntegration/reference/validation

valwfs=( "crt" "pds" "tpcsim" "tpcreco" )
nmiss=0

echo "looking for test references..."
for wf in ${valwfs[@]}; do

  if [ -e $refdir/test/$wf/ci_validation_histos.root ]; then
    echo -e "\033[0;32m \t$wf reference found \033[0m" #print this in green
  else
    echo -e "\033[0;31m \t$wf reference MISSING\033[0m" #print this in red
  fi

done

echo "looking for production references..."
for wf in ${valwfs[@]}; do

  if [ -e $refdir/$wf/ci_validation_histos.root ]; then
    echo -e "\033[0;32m \t$wf reference found \033[0m" #print this in green
  else
    echo -e "\033[0;31m \t$wf reference MISSING\033[0m" #print this in red
  fi

done

