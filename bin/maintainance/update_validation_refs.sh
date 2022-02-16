#!/bin/bash

# file: update_validation_refs.sh
# author: chilgenb@fnal.gov
# usage:  update_validation_refs.sh <CI build number>
# brief: given a CI build number, determine the validation workflow, code version,
#        next, check the persistent location for references and make changes as neccessary
#        to "retire" the existing reference file (unlink). then, move the new reference file
#        to the persistent location and create link to the new file
# args:
#   $1 -> (int) CI build number found on validation dashboard, usually 4 digits
#####################################################################################

### helper functions ###
file_exist(){
  if [ -f $1 ]; then
    echo -e "found reference file at \n\t $1"
  else
    echo "ERROR: could not locate reference file $1"
    kill -INT $$
  fi
}

# be careful not to overwrite files in output directory ($1)
check_outdir(){
  if [ `IFDH_FIND_MATCHING_DEPTH=1 ifdh findMatchingFiles ${1}/ ci_validation_histos.root | wc -l` == 0 ] && [ `IFDH_FIND_MATCHING_DEPTH=1 ifdh findMatchingFiles ${1}/ ci_validation_histos\*.root | wc -l` -gt 0 ]
  then
    echo "WARNING: validation files present in output directory, but no link (ci_validation_histos.root)."
    ifdh ls $1
    echo "CONTINUE? (y/n)"
    read usercont
    if [ $usercont != "y" ]; then
      echo "aborting..."
      kill -INT $$
    fi
  fi

  if [ `IFDH_FIND_MATCHING_DEPTH=1 ifdh findMatchingFiles ${1}/ ci_validation_histos_${tag}.root | wc -l` != 0 ]
  then
    echo "ERROR: a reference file for tag $tag already exists in the target directory! exiting..."
    kill -INT $$
  fi
}

# move reference file from scratch area ($1) to target directory ($2)
move() {

  if [ `IFDH_FIND_MATCHING_DEPTH=1 ifdh findMatchingFiles $2/ ci_validation_histos.root | wc -l` -gt 0 ] # is there a symlink there?
  then # dereference symlink and extract reference version

    fpath=`readlink -f $2/ci_validation_histos.root` #dereference symlink to get actual file path
    echo "reference path: $fpath"
    let beg=${#2}+1
    let size=${#fpath}-$beg
    fname=${fpath:$beg:$size}
    curtag=`cut -d 'v' -f3 <<< $fname`
    curtag=`cut -d '.' -f1 <<< $curtag`
    echo "found current reference file, $fname, with tag $curtag"

    echo "WARNING: your are about to retire the current ref histos in $2 with tag $curtag and replace them with $tag. do you wish to continue? (y/n)"
    read userin
    if [[ $userin != "y" ]]; then
      echo "  aborting..."
      kill -INT $$
    fi

    unlink $2/ci_validation_histos.root

  fi

  echo -e "\nmoving new references to $2"
  echo -e "\tifdh mv $1 ${2}/ci_validation_histos_${tag}.root"
  ifdh mv $1 ${2}/ci_validation_histos_${tag}.root
  here=`pwd`
  cd $2
  ln -s ci_validation_histos_${tag}.root ci_validation_histos.root # only relative symlinks work
  cd $here
}


####### main body #######
# check usage and environment
if [ ${#@} -eq 0 ]; then
  echo "usage: update_validation_refs.sh <CI build number>"
  exit 1
fi

if [ `ups active | grep ifdhc | wc -l` == 0 ]; then
  echo "ERROR: you must first do 'setup ifdhc'"
  exit 1
fi

# figure out if we're doing SBND or ICARUS stuff
source sbnci_setcodename.sh 
echo -e "\nlooking for ${expName}_ci/${1}..."

# look in scratch for a CI build with a matching number to the input and make sure it is unique
#   input dir structure is 
#      $indir/<validation WF>/<tag>/(CI_build_lar_ci_<build number> or just <build number>)/validation/
#   if subdir starts with CI_build_lar_ci it is a test submission
indir=/pnfs/${expName}/scratch/ci_validation/
nfound=`ifdh ls  $indir*/*/*$1/$1 recursion_depth==0 | wc -l`
if [ $nfound == 0 ]; then
  echo "ERROR: no build with number $1 was found"
  exit 1
elif [ $nfound -gt 1 ]; then
  echo "ERROR: multiple builds with number $1 were found"
  exit 1
fi

## determine validation workflow, code version, and if it's test or production
##   look for workflows and subdirectories (if applicable) defined below
valwfs=( "crt" "pds" "tpcsim" "tpcreco" )
tpcrecodirs=( "pfpvalidation" "pfpslicevalidation" "recoeff" "showervalidation" )
valwf=""

if [ `ifdh findMatchingFiles /pnfs/${expName}/scratch/ci_validation/*/*/*$1/validation *.root | wc -l` == 0 ]; then
  if [ `ifdh findMatchingFiles /pnfs/${expName}/scratch/ci_validation/*/*/*$1/$1/validation *.root | wc -l` == 0 ]; then
    echo "ERROR: could not locate CI build with number $1"
    exit 1

  else
    indir=`ifdh ls  $indir*/*/*$1/$1 recursion_depth==0`
  fi

else
  indir=`ifdh ls  $indir*/*/*$1 recursion_depth==0`
fi

echo -e "\nlooking for reference files in $indir"

for wf in ${valwfs[@]}; do
  if [[ $indir == *"$wf"* ]]; then
    echo "identified validation workflow, $wf"
    valwf="$wf"
    break
  fi
done

if [[ $valwf == "" ]]; then
  echo "ERROR: failed to identify validation workflow"
  kill -INT $$
fi

tag="$(cut -d '/' -f7 <<<$indir)"
echo "new reference tag is $tag"

# is it test or production?
outdir=/pnfs/${expName}/persistent/ContinuousIntegration/reference/validation
if [[ $indir == *"CI_build_lar_ci"*  ]]; then
  outdir="${outdir}/test"
fi
outdir=${outdir}/${valwf}

# done extracting info, move on to transfer/link update
# get absolute file path(s) to reference histogram(s) and move it(them) to 
# target directory in dCache persistent
refpathbase=${indir}validation/$valwf
if [[ $valwf == "tpcreco" ]]; then
  # loop over subdirectories
  for sub in ${tpcrecodirs[@]}; do
    refpath=$refpathbase/${sub}/ci_validation_histos.root
    file_exist $refpath
    check_outdir "${outdir}/${sub}"
    move $refpath "${outdir}/${sub}"
  done

# only tpcreco has subdirectories (for now)
else
  refpath=$refpathbase/ci_validation_histos.root
  file_exist $refpath
  check_outdir $outdir
  move $refpath $outdir
fi

echo "done. references updated successfully :)"
