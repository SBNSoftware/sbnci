#!/bin/bash

function get_proxy() {

  [ ${expName} ] || source sbnci_setcodename.sh

  echo "getting ${expName} analysis proxy"
  cigetcert -s 'fifebatch.fnal.gov'

  voms-proxy-info -exists -vo -valid 24:00 > /dev/null || 
      voms-proxy-init -noregen -rfc -voms "fermilab:/fermilab/${expName}/Role=Analysis" --valid 200:00

}

function warn_proxy() {
  echo "no valid proxy found! you must first run 'source get_proxy.sh'"
  exit 1
}

function check_proxy() {
  echo "checking proxy"
  voms-proxy-info || warn_proxy
}
