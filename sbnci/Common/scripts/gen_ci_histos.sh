#!/bin/bash

###################################################################
# flag options

a_flag='' # all WFs
p_flag='' # PDS sim/reco WF
r_flag='' # TPC reco WF
c_flag='' # CRT sim/reco WF
s_flag='' # TPC sim/cal WF

INFILE=''

print_usage() {
  printf "Usage: ./gen_ci_histos.sh <workflow options> -f <input file name>"
  printf "\n -a all validation workflows (overrides other options)"
  printf "\n -p PDS sim/reco validation workflow"
  printf "\n -r TPC reco validation workflow"
  printf "\n -c CRT sim/reco validation workflow"
  printf "\n -s TPC sim/cal validation workflow"
}

while getopts 'aprcsf:' flag; do
  case "${flag}" in
    a) a_flag='true' ;;
    p) p_flag='true' ;;
    r) r_flag='true' ;;
    c) c_flag='true' ;;
    s) s_flag='true' ;;
    f) INFILE="${OPTARG}" ;; 
    *) print_usage
       exit 1 ;;
  esac
done

#################################################################
# main

# Set code name and version to have experiment agnostic env vars
source sbnci_setcodename.sh

# execute plot scripts

if [[ $a_flag || $p_flag ]]
then
  plotScript="${SBNCI_DIR}/scripts/sbnciplot_pdsvalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)

elif [[ $a_flag || $r_flag ]]
then
  plotScript="${SBNCI_Dir}/scripts/sbnciplot_pfpvalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)
  
  plotScript="$SBNCI_DIR/scripts/sbnciplot_pfpslicevalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)
  
  plotScript="$SBNCI_DIR/scripts/sbnciplot_showervalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)
  
  plotScript="$SBNCI_DIR/scripts/sbnciplot_recoeff_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)

elif [[ $a_flag || $c_flag ]]
then
  plotScript="${SBNCI_DIR}/scripts/sbnciplot_crtvalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)

elif [[ $a_flag || $s_flag ]]
then
  plotScript="${SBNCI_DIR}/scripts/sbnciplot_tpcsimvalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)

fi
