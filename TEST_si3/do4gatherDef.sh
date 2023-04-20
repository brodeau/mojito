#!/bin/bash

#. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="${MOJITO_DIR}/gatherDef.py"

for RESKM in ${LCOARSEN}; do

    for NEMO_EXP in ${LIST_NEMO_EXP}; do
        echo; echo

        CMD="${EXE} ./npz nemoTsi3_idlSeed ${RESKM} NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}"

        echo
        echo " *** About to launch:"; echo "     ${CMD}"; echo

        ${CMD}


    done

done
