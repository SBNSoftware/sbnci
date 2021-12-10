#!/bin/bash

if [[ ${#1} == 0 || ${#2} != 0 ]]; then
  echo "usage: ./generate_bnb_intrnue_mix_validation_sample_test vXX_YY_ZZ"
  exit 1
fi 

source check_proxy.sh; check_proxy

source sbnci_setcodename.sh

if [[ $(ups list -aK+ $codeName $1 | wc -l) == 0 ]]; then
  echo "ERROR: could not find $codeName $1 in CVMFS!"
  exit 1
fi

echo "producing bnb_intrnue_mix validation (test run) sample for $codeName with tag $1"

trigger --build-delay 0 --jobname ${expName}_ci --workflow $expWF --gridwf-cfg cfg/${expName}/grid_workflow_${expName}_generate_bnb_intrnue_mix_validation_sample_test.cfg --revisions "SBNSoftware/$codeName@$1" --testmode
