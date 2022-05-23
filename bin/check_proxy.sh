#!/bin/bash

function get_proxy() {

  [ ${expName} ] || source sbnci_setcodename.sh

  echo "getting ${expName} analysis proxy"
  kx509 --minhours 24 #check cert has enough time left

  voms-proxy-info -exists -vo -valid 24:00 > /dev/null || 
      voms-proxy-init -noregen -rfc -voms "fermilab:/fermilab/${expName}/Role=Analysis" --valid 200:00

}

function warn_proxy() {
  echo "no valid proxy found! you must first run 'source get_proxy.sh'"
  exit 1
}

function check_proxy() {
  echo "checking proxy"
  voms-proxy-info -exists -vo -valid 24:00 || warn_proxy # check if we already have a proxy and it will last long enough
  cigetcert -s fifebatch.fnal.gov # check we have a valid proxy on myproxy server
}
