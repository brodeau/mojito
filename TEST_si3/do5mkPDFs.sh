#!/bin/bash

. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="${MOJITO_DIR}/mkpdfs.py"

DIR_NPZ="./npz"
#DIR_NPZ="./BAK_DEF"

for RESKM in ${LCOARSEN[*]}; do

    for NEMO_EXP in ${LIST_NEMO_EXP}; do

        echo; echo

        fdiv=`\ls ${DIR_NPZ}/def_DIV_NEMO-SI3_NANUK4_${NEMO_EXP}_${RESKM}km_${YEAR}????-${YEAR}????.npz`

        nbh=`echo ${fdiv} | wc -w`
        if [ ${nbh} -ne 1 ]; then
            if [ ${nbh} -eq 0 ]; then
                echo "ERROR: no divergence file found!!!"
            else
                echo "PROBLEM: more or less than 1 file found !!!"
            fi
            exit
        fi


        CMD="${EXE} ${fdiv}"
        echo "  ==> ${CMD}"; echo
        ${CMD}
        echo; echo


    done

done
