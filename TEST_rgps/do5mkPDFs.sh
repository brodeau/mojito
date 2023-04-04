#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/statistics.py"




fdiv=`\ls npz/def_DIV_RGPS_dt${DT_BINS_H}_${RESKM}km_${YEAR}????-${YEAR}????.npz`

nbh=`echo ${fdiv} | wc -w`
if [ ${nbh} -ne 1 ]; then
    echo "PROBLEM: more than 1 file found !!!"; exit
fi

CMD="${EXE} ${fdiv}"
echo "  ==> ${CMD}"; echo
${CMD}
echo; echo
