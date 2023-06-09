#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/gatherDef.py"

for NEMO_EXP in ${LIST_NEMO_EXP}; do
    echo; echo

    CMD="${EXE} ./npz ${DT_BINS_H} ${RESKM} NEMO-SI3_${NEMO_CONF}_${NEMO_EXP} ${USE_S}"

    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo

    ${CMD}


done
