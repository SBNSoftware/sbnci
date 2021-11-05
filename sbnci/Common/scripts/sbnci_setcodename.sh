#!/bin/bash

if [ ${SBNDCODE_VERSION} ]
then
  export codeName="sbndcode"
  export codeVersion=${SBNDCODE_VERSION}
  export expCIDir="${SBNCI_DIR}/SBND"
  export expName="sbnd"

elif [ ${ICARUSCODE_VERSION} ]
then
  export codeName="icaruscode"
  export codeVersion=${ICARUSCODE_VERSION}
  export expCIDir="${SBNCI_DIR}/ICARUS"
  export expName="icarus"

else
  echo "no active version of sbndcode or icaruscode found! (should never happen)"
fi
