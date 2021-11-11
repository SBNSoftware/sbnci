#!/bin/bash

if [ ${SBNDCODE_VERSION} ]
then
  export codeName="sbndcode"
  export codeVersion=${SBNDCODE_VERSION}
  export expCIDir="${SBNCI_DIR}/SBND"
  export expName="sbnd"
  export expWF="CI_VALIDATION_SBND"

elif [ ${ICARUSCODE_VERSION} ]
then
  export codeName="icaruscode"
  export codeVersion=${ICARUSCODE_VERSION}
  export expCIDir="${SBNCI_DIR}/ICARUS"
  export expName="icarus"
  export expWF="CI_VALIDATION_ICARUS"

else
  echo "no active version of sbndcode or icaruscode found! (should never happen)"
fi
