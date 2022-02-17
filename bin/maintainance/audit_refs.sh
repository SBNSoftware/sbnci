#!/bin/bash

# file: audit_refs.sh
# author: chilgenb@fnal.gov
# usage: audit_refs.sh
# brief: executable that checks that reference files are present for all validation WFs
# args: none 


# global vars
valwfs=( "crt" "pds" "tpcsim" "tpcreco" )
tpcrecodirs=( "pfpslicevalidation" "pfpvalidation" "recoeff" "showervalidation" )
nmiss=0
vers=0

# some useful functions
# given a path ($1) to a reference file link, extract the version
get_version(){
  fpath=`readlink -f $1` #dereference symlink to get actual file path
  #echo "reference path: $fpath"
  fname=`basename $fpath`
  
  curtag=`cut -d 'v' -f3 <<< $fname`
  curtag=`cut -d '.' -f1 <<< $curtag`
  if [ $curtag ]; then
    #curtag="?"
    curtag="v$curtag"
  else
    #curtag="v$curtag"
    curtag="?"
  fi

  if [ $vers != 0 ] && [ $vers != $curtag ]
  then
    echo -e "\033[0;31m \tWARNING: inconsistent reference versions:\033[0m $vers vs. $curtag"
  fi

  vers=$curtag
}


# check if a file exists and print a formatted result
# -- args --
# $1: full path to file
# $2: workflow name
# $3: subdirectory (for tpcreco only)
exist() {

  sub=""
  if [ $2 == "tpcreco" ]; then
    sub="/$3"
  fi

  if [ -e $1 ]; then
    get_version $1
    echo -e "\033[0;32m \t$2$sub reference found with version $vers \033[0m" #print this in green
  else
    echo -e "\033[0;31m \t$2$sub reference MISSING\033[0m" #print this in red
    let nmiss=$nmiss+1
  fi

}

# look for ref files for all validation WFs base directory
# $1: full directory path
look() {

  for wf in ${valwfs[@]}; do

    if [ $wf == "tpcreco" ]; then

      for sub in ${tpcrecodirs[@]}; do

        path=$1/$wf/$sub/ci_validation_histos.root
        exist $path $wf $sub
      done

    else
      path=$1/$wf/ci_validation_histos.root
      exist $path $wf
    fi
  
  done
}

############################################################################
# main body

source sbnci_setcodename.sh # figure out if we're doing SBND or ICARUS stuff
refdir=/pnfs/${expName}/persistent/ContinuousIntegration/reference/validation

echo -e "\nlooking for files in ${refdir}..."

echo -e "\nchecking test references..."
look "$refdir/test"

echo -e "\nchecking production references..."
look "$refdir"

if [ $nmiss == 0 ]; then
  echo -e "\033[0;32m \nall references found :) \033[0m" #print this in green
elif [ $nmiss == 1 ]; then
  echo -e "\033[0;31m \n1 reference MISSING\033[0m" #print this in red
else
  echo -e "\033[0;31m \n$nmiss references MISSING\033[0m" #print this in red
fi

echo -e "\naudit complete."
