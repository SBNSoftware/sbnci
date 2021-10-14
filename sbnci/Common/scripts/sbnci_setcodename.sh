#!/bin/bash

if [ ${SBNDCODE_VERSION} ]
then
  export codeName="sbndcode"
  export codeVersion=${SBNDCODE_VERSION}
  export expCIDir="${SBNCI_DIR}/SBND"

elif [ ${ICARUSCODE_VERSION} ]
then
  export codeName="icaruscode"
  export codeVersion=${ICARUSCODE_VERSION}
  export expCIDir="${SBNCI_DIR}/ICARUS"

else
  echo "no active version of sbndcode or icaruscode found! (should never happen)"
fi
