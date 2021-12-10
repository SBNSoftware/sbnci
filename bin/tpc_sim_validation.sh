#!/bin/bash

if [ ${#@} == 0 ]; then
  echo "usage: tpc_sim_validation.sh space-separated list of repo@branch's (SBNSoftware repos only)"
  exit 1
fi 

source check_proxy.sh; check_proxy

echo "preparing to validate ${#@} branches"

# N.B. this only works if all repos are in SBNSoftware
branchstr=""
for branch in "$@"
do
  branchstr="SBNSoftware/$branch $branchstr"
done

source sbnci_setcodename.sh

echo "validating TPC sim for $branchstr"

trigger --build-delay 0 --jobname ${expName}_ci --workflow $expWF --gridwf-cfg cfg/${expName}/grid_workflow_${expName}_tpcsim_all.cfg --revisions "$branchstr"
