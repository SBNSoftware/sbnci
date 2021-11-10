#!/bin/bash

function get_proxy() {

  echo "getting proxy"
  setup cigetcert
  cigetcert -s 'fifebatch.fnal.gov'
  voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'
  voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/icarus/Role=Analysis'

}

function warn_proxy() {
  echo "no valid proxy found! you must first run 'source get_proxy.sh'"
  exit 1
}

function check_proxy() {
  echo "checking proxy"
  voms-proxy-info || warn_proxy
}
