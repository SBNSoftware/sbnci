#!/bin/bash

# file: check_branches.sh
# created: 16 Dec 2021
# author: Chris Hilgenberg (chilgenb@fnal.gov)
# brief: this script takes a string of revisions to be passed
#        to the CI system and checks that the format is correct
if [ ${#1} == 0 ]; then
  echo "usage: check_branches.sh space-separated list of repo@branch's (SBNSoftware repos only)"
  exit 1
fi

# branches should have the form 'SBNSoftware/<repo>@<tag or branch>'
#   e.g. SBNSoftware/sbndcode@vXX_YY_ZZ
#   if the branch is not valid, CI will take develop branch instead
#   of either sbndcode or icaruscode
base='SBNSoftware/'
for branch in "$@"; do

  if [[ $branch != *"@"* ]]; then
    echo "ERROR: incorrect branch syntax"
    echo -e "\toffending branch: $branch"
    exit 1
  fi

  # trim off 'SBNSoftware' to get repo@tag
  if [[ $branch == "$base"* ]]; then
    branch=${branch:${#base}:${#branch}}
  fi

  # repository name should not have '/'s
  repo=`echo $branch | cut -d'@' -f 1`
  if [[ "$repo" == *"/"* ]]; then
    echo "ERROR: invalid repository name $repo"
    exit 1
  fi
  
done
