#!/bin/bash

###################################################################
# flag options

all_flag='' # all WFs
pds_flag='' # PDS sim/reco WF

tpcreco_flag='' # TPC reco WF
pfp_flag=''
pfpslice_flag=''
recoeff_flag=''
shower_flag=''

crt_flag='' # CRT sim/reco WF
tpcsim_flag='' # TPC sim/cal WF


INFILE=''

print_usage() {
  printf "Usage: ./gen_ci_histos.sh <workflow options> -f <input file name>"
  printf "\n -a, --all         all validation workflows"
  printf "\n -p, --pds         PDS sim/reco validation workflow"
  printf "\n -r, --tpcreco     TPC reco validation workflow"
  printf "\n --pfp             PFP validation workflow"
  printf "\n --pfpslice        PFP slice validation workflow"
  printf "\n --recoeff         TPC reco efficiency validation workflow"
  printf "\n --shower          Shower reco efficiency validation workflow"
  printf "\n -c, --crt         CRT sim/reco validation workflow"
  printf "\n -s, --tpcsim      TPC sim/cal validation workflow"
  printf "\n"
}

while test $# -gt 0; do
  case "$1" in
    -a|--all)
      shift
      all_flag='true'
      ;;
    -c|--crt)
      shift
      crt_flag='true'
      ;;
    -p|--pds)
      shift
      pds_flag='true'
      ;;
    -s|--tpcsim)
      shift
      tpcsim_flag='true'
      ;;
    -r|--tpcreco)
      shift
      tpcreco_flag='true'
      ;;
    --pfp)
      shift
      pfp_flag='true'
      ;;
    --pfpslice)
      shift
      pfpslice_flag='true'
      ;;
    --recoeff)
      shift
      recoeff_flag='true'
      ;;
    --shower)
      shift
      shower_flag='true'
      ;; 
    -f|--file)
      INFILE="$2"
      shift 2
      ;;
    *) print_usage
       exit 1
       ;;
  esac
done

if [ $INFILE == '' ]
then
  echo -e "no input file provided (required)\n"
  print_usage
  exit 1
fi

#################################################################
# main

# Set code name and version to have experiment agnostic env vars
source sbnci_setcodename.sh

# execute plot scripts

if [[ $all_flag || $pds_flag ]]
then
  plotScript="${SBNCI_DIR}/scripts/sbnciplot_pdsvalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)
fi

if [[ $all_flag || $tpcreco_flag ]]
then
  if [[ $tpcreco_flag || $pfp_flag ]]
  then
    plotScript="${SBNCI_Dir}/scripts/sbnciplot_pfpvalidation_${expName}.C"
    root -l -b -q $plotScript\(\"${INFILE}\"\)
  fi

  if [[ $tpcreco_flag || $pfpslice_flag ]]
  then
    plotScript="$SBNCI_DIR/scripts/sbnciplot_pfpslicevalidation_${expName}.C"
    root -l -b -q $plotScript\(\"${INFILE}\"\)
  fi

  if [[ $tpcreco_flag || $shower_flag ]]
  then
    plotScript="$SBNCI_DIR/scripts/sbnciplot_showervalidation_${expName}.C"
    root -l -b -q $plotScript\(\"${INFILE}\"\)
  fi

  if [[ $tpcreco_flag || $recoeff_flag ]]
  then
    plotScript="$SBNCI_DIR/scripts/sbnciplot_recoeff_${expName}.C"
    root -l -b -q $plotScript\(\"${INFILE}\"\)
  fi
fi

if [[ $all_flag || $crt_flag ]]
then
  plotScript="${SBNCI_DIR}/scripts/sbnciplot_crtvalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)
fi

if [[ $all_flag || $tpcsim_flag ]]
then
  plotScript="${SBNCI_DIR}/scripts/sbnciplot_tpcsimvalidation_${expName}.C"
  root -l -b -q $plotScript\(\"${INFILE}\"\)
fi

