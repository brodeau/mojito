#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/diags/mkpdfs.py"

fdiv=`\ls npz/def_DIV_RGPS_dt${DT_BINS_H}_${RESKM}km.npz`

nbh=`echo ${fdiv} | wc -w`
if [ ${nbh} -ne 1 ]; then
    if [ ${nbh} -eq 0 ]; then
        echo "ERROR: no divergence file found!!!"
    else
        echo "PROBLEM: more or less than 1 file found !!!"
    fi
    exit
fi


CMD="${EXE} ${fdiv} ${MODE}"
echo "  ==> ${CMD}"; echo
${CMD}
echo; echo
