#!/bin/bash

#. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="${MOJITO_DIR}/gatherDef.py"

#DIR_NPZ="./npz"
DIR_NPZ="./BAK_DEF"


for RESKM in ${LCOARSEN}; do

    for NEMO_EXP in ${LIST_NEMO_EXP}; do
        echo; echo

        CMD="${EXE} ${DIR_NPZ}/${RESKM} nemoTsi3_NoBin ${RESKM} NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}"

        echo
        echo " *** About to launch:"; echo "     ${CMD}"; echo

        ${CMD}


    done

done
