#!/bin/bash

# @brief: wrapper script for lar_ci trigger command 
# @usage: validate <workflow> --ref <vXX_YY_ZZ or current> --revisions <repo@branches> --test (optional - use lar_ci "test mode")

# global vars
versc=""
versr=""
revs=""
branchstr=""
testmode=""
SBNCI_REF_VERSION="" # this variable needs to be exported to lar_ci (w/matching cfg definition)
workflow=""

GetHelp() {
  echo -e "\n validate.sh is a convenience script to make submitting validation jobs easier"
  echo -e "\n\tusage: validate <workflow> --ref <vXX_YY_ZZ or current> --revisions <repo@branches> --test (optional - use lar_ci 'test mode')"
  echo -e "\n\trequired parameters:"
  echo -e "\t\tworkflow is one of 'crt' , 'tpcreco', 'tpcsim' or 'pds'"
  echo -e "\t\t-r or --ref is vXX_YY_ZZ or use -c or --current (derives from ${expName}code@develop), the reference version of the code used for the validation"
  echo -e "\t\t--revisions is space-seperated list of SBNSoftware branches being tested in the form repo0@branch0 repo1@branch1 ..."
  echo -e "\t\t\tfor example --revisions ${expName}code@my_${expName}code_branch ${expName}util@vUU_VV_YY"
  echo -e "\n\toptional parameters:"
  echo -e "\t\t-t or --testmode: if specified, use lar_ci testmode (should be used before full validation)"
} 

CheckValidationWorkflow() {
  case "$1" in
    "crt")
      workflow="crt" ;;
    "tpcreco")
      workflow="mc_reco_all" ;;
    "tpcsim")
      workflow="tpcsim" ;;
    "pds")
      workflow="pds" ;;
  esac

  if [ "$workflow" == "" ]; then
    echo "ERROR: validation workflow '$1' is not recognized. Choose from the following options."
    echo "  'crt' , 'tpcreco', 'tpcsim' or 'pds'"
  fi
}

# helper functions
CompleteSbnsoftName(){
  for branch in $revs
  do
    if [ "$branchstr" == "" ]; then
      branchstr="\"SBNSoftware/$branch"
    else
      branchstr="SBNSoftware/$branch $branchstr"
    fi
  done
  branchstr="${branchstr}\""
  echo -e "\ncomplete branch names for list of revisions: ${branchstr}"
}

CheckReferenceVersion() {

  # must use grep -w to match to entier string delimited by whitespace
  nmatch=`cat /pnfs/${expName}/persistent/ContinuousIntegration/approved_reference_versions.txt | grep -w $1 | wc -l`

  if [ $nmatch == 0 ]; then
    echo "ERROR: reference version specified by user ($1) is not approved."
    echo "Choose from (# reference version # reference alias #)"
    cat /pnfs/sbnd/persistent/ContinuousIntegration/approved_reference_versions.txt

  elif [ $nmatch -gt 1 ]; then
    echo "ERROR: multiple matches ($nmatch) found for reference version $1"

  else
    str=$(cat /pnfs/${expName}/persistent/ContinuousIntegration/approved_reference_versions.txt | grep -w $1)
    strarr=($str)
    if [ "$1" == "${strarr[0]}" ]; then # $1 is a version
      if [ "${strarr[1]}" == "current" ]; then
        SBNCI_REF_VERSION="current"
      else
        SBNCI_REF_VERSION="$1"
      fi

    elif [ "$1" == "${strarr[1]}" ]; then # $1 is a version alias
      if [ "$1" == "current" ]; then
        SBNCI_REF_VERSION="current"
      else
        SBNCI_REF_VERSION="${strarr[0]}"
      fi

    else # shouldn't happen ever, but hey doesn't hurt to check
      echo "ERROR: $1 found in file but no match to version or alias fields"
    fi
  fi
}

# parse input args
SetReferenceArgs(){

  if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    GetHelp

  elif [ "$1" == "" ]; then
    echo "usage: validate.sh <workflow> --ref <vXX_YY_ZZ or current> --revisions <repo@branches> --test (optional)"
    echo "try 'validate.sh -h' for more information"

  else
    #echo "call SetReferenceArgs with input $@" # for debugging
    CheckValidationWorkflow "$1"

    while [ "${1:-}" != "" ]; do
      case "$1" in
        "-c" | "--current")
          versc="current" ;;
        "-r" | "--ref")
          shift
          versr=$1 ;;
        "-t" | "--test")
          testmode="--testmode";;
        "--revisions")
          shift
          while [ "$1" != "" ]; do
            if [ "$revs" == "" ]; then
              revs="$1"
            else
              revs="$1 $revs"
            fi
 
            if [ "${1:0:1}" == "-" ] || [ "${1:0:2}" != "--" ]; then
              break
            else
              shift
            fi

          done
          #echo "revisions: '$revs'" # for debugging
          CompleteSbnsoftName "$revs" ;;
      esac
      shift
    done

    # check inputs
    if [ "$versc" != "" ]  && [ "$versr" != "" ]; then
      echo "ERROR: You cannot specify both --ref and --current"

    elif [ "$versc" == "" ] && [ "$versr" == "" ]; then
      echo "ERROR: You must specify a reference version: --ref vXX_YY_ZZ or --current"

    elif [ "$branchstr" == "" ]; then
      echo "ERROR: No revisions provided. You must provide --revisions <repo0@branch0 repo1@branch1 ...>"

    elif [ "$workflow" != "" ]; then

      # hack to handle inconsistent workflow naming 
      if [ "$testmode" == "" ]; then #testmode disabled
        if [ "$workflow" != "mc_reco_all" ]; then
          workflow="${workflow}_all"
        fi
      else #testmode enabled
        workflow="${workflow}_test"
      fi

      if [ "$versc" != "" ]; then
        SBNCI_REF_VERSION=$versc
      elif [ "$versr" != "" ]; then
        SBNCI_REF_VERSION=$versr
      fi

      CheckReferenceVersion $SBNCI_REF_VERSION
      echo "using reference version $SBNCI_REF_VERSION"

    fi # workflow is ok
  fi # not asking for help
}

# main body

 source sbnci_setcodename.sh # for releases (UNCOMMENT WHEN COMMITTING!)
 #source $MRB_INSTALL/sbnci/$MRB_PROJECT_VERSION/slf7.x86_64.e20.prof/bin/sbnci_setcodename.sh # for local development (COMMENT OUT WHEN COMMITTING!)
 SetReferenceArgs $@

 if [ "$SBNCI_REF_VERSION" != "" ]; then

   # hack to avoid sbndcode CI test failure with CAF checks with SBN2022A
   citests=""
   if [ "$SBNCI_REF_VERSION" == "v09_37_02_01" ] && [ "${expName}" == "sbnd" ]; then
     citests="--ci-tests \"nucosmics_g4_quick_test_sbndcode single_g4_quick_test_sbndcode single_reco2_quick_test_sbndcode compilation_test_sbndcode nucosmics_detsim_quick_test_sbndcode nucosmics_g4_quick_test_sbndcode nucosmics_gen_quick_test_sbndcode nucosmics_reco1_quick_test_sbndcode nucosmics_reco2_quick_test_sbndcode single_gen_quick_test_sbndcode single_reco1_quick_test_sbndcode\""
   fi

   gridwf="cfg/${expName}"
   if [ "$SBNCI_REF_VERSION" == "current" ]; then
     gridwf="${gridwf}/current"
   fi
   gridwf="${gridwf}/grid_workflow_${expName}_${workflow}.cfg"

   extras=""
   if [ "${citests}" != "" ]; then
     extras="${citests}"
   fi
   if [ "${testmode}" != "" ]; then
     if [ "$extras" != "" ]; then
       extras="${extras} ${testmode}"
     else
       extras="${testmode}"
     fi
   fi
   echo "trigger --build-delay 0 --jobname ${expName}_ci --workflow $expWF --gridwf-cfg $gridwf --revisions $branchstr -e SBNCI_REF_VERSION=$SBNCI_REF_VERSION $extras"

   cmd="trigger --build-delay 0 --jobname ${expName}_ci --workflow $expWF --gridwf-cfg $gridwf --revisions $branchstr -e SBNCI_REF_VERSION=$SBNCI_REF_VERSION $extras"
   eval $cmd

   if [ "$testmode" != "" ]; then 
     echo -e "\ntest validation jobs submitted. go to the test dashboard to view your results."
     echo -e "\thttps://dbweb9.fnal.gov:8443/TestCI/app/ns:${expName}/view_builds\n"

   else
     echo -e "\nvalidation jobs submitted. go to the dashboard to view your results."
     echo -e "\thttps://dbweb9.fnal.gov:8443/LarCI/app/ns:${expName}/view_builds\n"
   fi

 else
   echo "ERROR: no viable reference version could be deduced"
 fi



