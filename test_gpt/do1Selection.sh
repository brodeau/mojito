#!/bin/bash

. ./conf.bash

EXE="python3 -u ${MOJITO_DIR}/rgps/01_selection_xy.py"

YEAR=`echo ${DATE1} | cut -c1-4`

#if [ `echo ${DATE2} | cut -c1-4` -ne ${YEAR} ]; then
#    echo "ERROR: DATE1 and DATE2 have 2 different years!" ; exit
#fi

MMDD1=`echo ${DATE1} | cut -c5-8`
MMDD2=`echo ${DATE2} | cut -c5-8`

if [ ${RESKM} -ge 50 ] && [ ${DT_BINS_H} -ge 72 ]; then
    CMD="${EXE} ${RGPS_DIR}/${FILIN} ${DATE1} ${DATE2} ${DT_BINS_H} $((DT_BINS_H*2/3))"
else
    CMD="${EXE} ${RGPS_DIR}/${FILIN} ${DATE1} ${DATE2} ${DT_BINS_H}"
fi

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
