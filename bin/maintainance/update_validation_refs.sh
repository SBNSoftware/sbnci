#!/bin/bash

if [ ${#@} -eq 0 ]; then
  echo "usage: ./update_validation_refs.sh <CI build number>"
  exit 1
fi

if [ `ups active | grep ifdhc | wc -l` == 0 ]; then
  echo "ERROR: you must first do 'setup ifdhc'"
  exit 1
fi

source sbnci_setcodename.sh # figure out if we're doing SBND or ICARUS stuff
echo "looking for ${expName}_ci/$1"


# input dir structure is $indir/<validation WF>/<tag>/(CI_build_lar_ci_<build number> or just <build number>)/validation/
#  if subdir starts with CI_build_lar_ci it is a test submission
indir=/pnfs/${expName}/scratch/ci_validation/
nfound=`ifdh ls  $indir*/*/*$1 recursion_depth==0 | wc -l`
if [ $nfound == 0 ]; then
  echo "ERROR: no build with number $1 was found"
  exit 1
elif [ $nfound -gt 1 ]; then
  echo "ERROR: multiple builds with number $1 were found"
  exit 1
fi

valwfs=( "crt" "pds" "tpcsim" "tpcreco" )
tpcrecowfs=( "pfpvalidation" "pfpslicevalidation" "recoeff" "showervalidation" )

indir=`ifdh ls  $indir*/*/*$1 recursion_depth==0`
echo "looking for reference files in $indir"

valwf=""
for wf in ${valwfs[@]}; do
  if [[ $indir == *"$wf"* ]]; then
    echo "identified validation workflow, $wf"
    valwf="$wf"
    break
  fi
done

if [[ $valwf == "" ]]; then
  echo "ERROR: failed to identify validation workflow"
  exit 1
fi

tag="$(cut -d '/' -f7 <<<$indir)"
echo "new reference tag is $tag"

# get absolute file path(s) to reference histogram(s)
refpathbase=${indir}validation/$valwf
if [[ $valwf == "tpcreco" ]]; then
  for wf in ${tpcrecowfs[@]}; do
    refpath=$refpathbase/$wf/ci_validation_histos.root

    if [ -f $refpath ]; then
      echo -e "found reference file at \n\t $refpath"
    else
      echo "ERROR: could not locate reference file"
      exit 1
    fi

  done

else
  refpath=$refpathbase/ci_validation_histos.root

  if [ -f $refpath ]; then
    echo -e "found reference file at \n\t $refpath"
  else
    echo "ERROR: could not locate reference file"
    exit 1
  fi
fi

# is it test or production?
outdir=/pnfs/${expName}/persistent/ContinuousIntegration/reference/validation
if [[ $indir == *"CI_build_lar_ci"*  ]]; then
  outdir="${outdir}/test"
fi

# be careful not to overwrite files
outdir=$outdir/$valwf
if [ `ifdh ls $outdir/ci_validation_histos.root | wc -l` -gt 0 ]; then
  echo "WARNING: your are about to retire the current ref histos. do you wish to continue? (y/n)"
  read userin
  if [[ $userin != "y" ]]; then
    echo "  aborting..."
    exit 1
  fi

  unlink $outdir/ci_validation_histos.root

else
  echo "moving new refs to $outdir"
  ifdh mv $refpath ${outdir}/ci_validation_histos_${tag}.root
  ln -s ${outdir}/ci_validation_histos_${tag}.root ${outdir}/ci_validation_histos.root
fi

echo "references updated successfully"

exit 0
