#!/bin/bash

# @brief: wrapper script for lar_ci trigger command 
# @usage: validate <workflow> --ref <vXX_YY_ZZ or current> --revisions <repo@branches> --test (optional - use lar_ci "test mode")

# global vars
versc=""
versr=""
revs=""
branchstr=""
extrabranchstr=""
testmode=""
SBNCI_REF_VERSION="" # this variable needs to be exported to lar_ci (w/matching cfg definition)
workflow=""
testtag=""
reftag=""
currenttag=""
gridwf=""

# three possible tag types (Chris' probably wrong naming)
t1="vXX_YY_ZZ" # version core
t2="${t1}_ii"  # hotfix
t3="${t2}pjj" # hotfix patch


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

  repolist="/pnfs/${expName}/persistent/ContinuousIntegration/supported_repos.txt"

  for br in ${revs[@]}
  do

    # parse org/repo@branch elements
    repo=$(echo $br | cut -d '@' -f 1)
    branch=$(echo $br | cut -d '@' -f 2)
    org=$(echo $repo | cut -s -d '/' -f 1)
    base=$repo
    if [ "$org" != "" ]; then
      base=$(echo $repo | cut -d '/' -f 2)
    fi

    if [ "$base" == "${expName}code" ] || [ "$base" == "${expName}util" ] ; then
      if [ "$branchstr" == "" ]; then
        branchstr="SBNSoftware/$base@$branch"
      else
        branchstr="$branchstr SBNSoftware/$base@$branch"
      fi
    fi

    repoinfile=$(cat $repolist | grep $repo) #gets line from repolist in format org/repo
    if [ "$repoinfile" == "" ]; then
      echo "ERROR: branch $br is not recognized. Contact the SBN validation team if you need $repo added."
      return
    fi

    larste="LARSOFT_SUITE_v"
    # if it's a larsoft repo, it's either in LArSoft or a fork in SBNSoftware
    if [ "${base:0:3}" == "lar" ]; then 

      if [ "${branch:0:$(expr ${#larste})}" == "$larste" ]; then #if larsuite in branch name

        if [ "$org" == "SBNSoftware" ]; then  # this tag doesn't exist there
          echo "ERROR: LArSoft-exclusive tag mixed with SBNSoftware"
          return
        fi

        org="LArSoft"

      elif [ "$(cat $repolist | grep "SBNSoftware/$base")" != "" ]; then

        org="SBNSoftware"

      else

        echo "ERROR: branch $branch from LArSoft/ should be of the form LARSOFT_SUITE_vAB_CD_EF"
        return

      fi

    else #if not a larsoft repo, it's an SBN repo
      org="SBNSoftware"

    fi


    if [ "$base" != "${expName}code" ] || [ "$base" != "${expName}util" ] ; then
      if [ "$extrabranchstr" == "" ]; then
        extrabranchstr="$org/$base@$branch"
      else
        extrabranchstr="$extrabranchstr $org/$base@$branch"
      fi
    fi   

  done

  branchstr="\"$branchstr\"" #add quotes
  extrabranchstr="\"$extrabranchstr\""

  #echo -e "\ncomplete branch names for list of revisions: ${branchstr}"
  #echo -e "\ncomplete branch names for list of extra revisions: ${extrabranchstr}"

}

#  In order to determine the tag associated with the branch(es) being tested,
#  we need to create a work area, clone the <exp>code repo, checkout the branch
#  of that repo being tested, then use git commands to extract the tag.
GetTestTag() {

  # $revs does not have SBNSoftware in branch name
  if [ "$revs" == "" ]; then
    echo "ERROR: GetTestTag cannot extract tag from empty branch list"

  else
    expbranch=""
    brancharr=( $revs )
    for branch in ${brancharr[@]}; do
      if [ "$(echo "$branch" | grep "${expName}code")" != "" ]; then
        expbranch="$branch"
        break
      fi
    done
    expbranch=`echo $expbranch | cut -d '/' -f 2` # if org in branch name, cut it out

    if [ "$expbranch" == "" ]; then
      echo "ERROR: GetTestTag could not identify the experiment code branch."
      echo "WARNING: your branch(es) might not be tested with the appropriate sim/reco workflow"

    else
      exprepo="${expName}code"
      expbranch="${expbranch:$(expr ${#exprepo}+1):${#expbranch}}"

      tmpdir=/${expName}/app/users/ContinuousIntegration/scratch/`date +"%s"`
      if [ -e $tmpdir ]; then
        echo "ERROR: GetTestTag failed to create test directory"
      else
        here=`pwd`
        mkdir $tmpdir
        cd $tmpdir

        git clone --branch "$expbranch" "https://github.com/SBNSoftware/${exprepo}.git" 1> /dev/null 2> /dev/null
        if [ -d "$exprepo" ]; then
          cd ${exprepo}
        else
          echo "ERROR: git clone failed! unable to deduce test version"
          return
        fi

        testtag=`git describe --tags --abbrev=0`

        cd $here
        rm -rf $tmpdir
      fi
    fi
  fi
}

# pass a tag ($1) to make sure it's sane and get it's value as a number
CheckTagFormat(){
  if [ "${1:0:1}" != "v" ]; then
    echo "WARNING: strange tag ($1 doesn't start with 'v') detected! test sim/reco workflow may not be correct"
  fi

  if [ ${#1} != ${#t1} ] && [ ${#1} != ${#t2} ] && [ ${#1} != ${#t3} ]; then
    echo "WARNING: strange tag ($1 is not of the form vXX_YY_ZZ_iipjj) detected! test sim/reco workflow may not be correct"
  fi

  
}

# read from text file a list of approved reference tags and make sure
# that reference tag passed by user is on the list
CheckReferenceVersion() {

  # must use grep -w to match to entier string delimited by whitespace
  nmatch=`cat /pnfs/${expName}/persistent/ContinuousIntegration/approved_reference_versions.txt | grep -w $1 | wc -l`

  if [ $nmatch == 0 ]; then
    echo "ERROR: reference version specified by user ($1) is not approved."
    echo "Choose from (# reference version # reference alias #)"
    cat /pnfs/${expName}/persistent/ContinuousIntegration/approved_reference_versions.txt

  elif [ $nmatch -gt 1 ]; then
    echo "ERROR: multiple matches ($nmatch) found for reference version $1"

  else
    str=$(cat /pnfs/${expName}/persistent/ContinuousIntegration/approved_reference_versions.txt | grep -w $1)
    strarr=($str)
    if [ "$1" == "${strarr[0]}" ]; then # $1 is a version

      reftag="$1"

      if [ "${strarr[1]}" == "current" ]; then
        SBNCI_REF_VERSION="current"
      else
        SBNCI_REF_VERSION="$1"
      fi

    elif [ "$1" == "${strarr[1]}" ]; then # $1 is a version alias

      reftag="${strarr[0]}"

      if [ "$1" == "current" ]; then
        SBNCI_REF_VERSION="current"
      else
        SBNCI_REF_VERSION="${strarr[0]}"
      fi

    else # shouldn't happen ever, but hey doesn't hurt to check
      echo "ERROR: $1 found in file but no match to version or alias fields"
    fi

    # also need to know what tag 'current' points to
    str=$(cat /pnfs/${expName}/persistent/ContinuousIntegration/approved_reference_versions.txt | grep -w current)
    strarr=($str)
    currenttag="${strarr[0]}"

  fi

  # decide what sim/reco workflow to use by comparing the tag the test branch 
  #  is based off of to that of the specified reference as well as that of 'current'
  # TODO support arbitrary number of sim/reco workflows instead of 
  #  just current vs. most recent prod branch - need to add more subdirs to
  #  lar_ci/cfg/<exp>/ similarly to how was done for 'current'
  #echo "current tag: $currenttag"
  CheckTagFormat "$currenttag"
  CheckTagFormat "$reftag"
  GetTestTag
  CheckTagFormat "$testtag"

  # convert tags to six-digit integers assuming hot fixes and hot fix-patches
  #  don't confuse things (should only have one instance per production branch/current)
  expr testnum="${testtag:1:2}${testtag:4:2}${testtag:7:2}" > /dev/null
  expr refnum="${reftag:1:2}${reftag:4:2}${reftag:7:2}" > /dev/null
  expr curnum="${currenttag:1:2}${currenttag:4:2}${currenttag:7:2}" > /dev/null

  gridwf="cfg/${expName}"

  if [ $testnum -gt $refnum ]; then
    gridwf="${gridwf}/current"

  # use this for more than two options (see TODO above)
  #elif [ "$testnum" -lt "$curnum" ] && [ "$testnum" -ge "$refnum" ]; then

  elif [ "${testnum}" -lt "$refnum" ]; then
    echo "ERROR: CI system does not support validating branches based on tags older than the reference version"
    echo "WARNING: using sim/reco workflow for 'current' which may not be compatible with your branch"
    gridwf="${gridwf}/current"

  fi

  gridwf="${gridwf}/grid_workflow_${expName}_${workflow}.cfg"
  
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
        "--revisions")
          shift
          while [ "$1" != "" ]; do
            if [ "$revs" == "" ]; then
              revs="$1"
            else
              revs="$1 $revs"
            fi
            if [ "${2:0:1}" == "-" ]; then
              break
            fi

            shift
          done
          #echo "revisions: '$revs'" # for debugging
          CompleteSbnsoftName "$revs" ;;

        "-c" | "--current")
          versc="current" ;;
        "-r" | "--ref")
          shift
          versr=$1 ;;
        "-t" | "--test")
          testmode="--testmode";;
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

      #echo "passing SBNCI_REF_VERSION=$SBNCI_REF_VERSION to CheckReferenceVersion"
      CheckReferenceVersion $SBNCI_REF_VERSION
      #echo "using reference version $SBNCI_REF_VERSION"

    fi # workflow is ok
  fi # not asking for help
}


# main body

 source sbnci_setcodename.sh # for releases (UNCOMMENT WHEN COMMITTING!)
 #source $MRB_INSTALL/sbnci/$MRB_PROJECT_VERSION/slf7.x86_64.e20.prof/bin/sbnci_setcodename.sh # for local development (COMMENT OUT WHEN COMMITTING!)
 SetReferenceArgs $@

 # TODO check that input branch list does not contain two branches from the same repo
 # set the correct environment var
 envextra=""
 if [ "$expName" == "sbnd" ]; then
   envextra="-e SBNDmodules_extra"

 elif [ "$expName" == "icarus" ]; then
   envextra="-e ICARUSmodules_extra"

 fi

 if [ "$SBNCI_REF_VERSION" != "" ]; then

   # hack to avoid sbndcode CI test failure with CAF checks with SBN2022A
   citests=""
   if [ "$SBNCI_REF_VERSION" == "v09_37_02_01" ] && [ "${expName}" == "sbnd" ]; then
     citests="--ci-tests \"nucosmics_g4_quick_test_sbndcode single_g4_quick_test_sbndcode single_reco2_quick_test_sbndcode compilation_test_sbndcode nucosmics_detsim_quick_test_sbndcode nucosmics_g4_quick_test_sbndcode nucosmics_gen_quick_test_sbndcode nucosmics_reco1_quick_test_sbndcode nucosmics_reco2_quick_test_sbndcode single_gen_quick_test_sbndcode single_reco1_quick_test_sbndcode\""
   fi

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
   if [ "$extrabranchstr" != "" ]; then
     if [ "$extras" != "" ]; then
       extras="${extras} ${envextra}=${extrabranchstr}"
     else
       extras="${envextra}=${extrabranchstr}"
     fi

   fi

   cmd="trigger --build-delay 0 --jobname ${expName}_ci --workflow $expWF --gridwf-cfg $gridwf --revisions $branchstr -e SBNCI_REF_VERSION=$SBNCI_REF_VERSION $extras"

   #echo "$cmd"

   eval $cmd 1>/dev/null #2>/dev/null

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



