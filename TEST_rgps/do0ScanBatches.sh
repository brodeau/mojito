#!/bin/bash

. ./conf.bash

EXE1="python3 -u ${MOJITO_DIR}/rgps/00_scan_batches.py"
EXE2="python3 -u ${MOJITO_DIR}/rgps/01_selection_xy.py"

YEAR=`echo ${DATE1} | cut -c1-4`

#if [ `echo ${DATE2} | cut -c1-4` -ne ${YEAR} ]; then
#    echo "ERROR: DATE1 and DATE2 have 2 different years!" ; exit
#fi

MMDD1=`echo ${DATE1} | cut -c5-8`
MMDD2=`echo ${DATE2} | cut -c5-8`



fout="./npz/RGPS_batch_selection_dt${DT_BINS_H}h_${DATE1}_${DATE2}_mode-${MODE}.npz"

if [ ! -f ${fout} ]; then

    CMD="${EXE1} ${RGPS_DIR}/${FILIN} ${DATE1} ${DATE2} ${DT_BINS_H} ${MODE}"
    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo
    ${CMD}
    echo; echo
    
fi

# Second script:
if [ ${RESKM} -ge 50 ] && [ ${DT_BINS_H} -ge 72 ]; then
    CMD="${EXE2} ${RGPS_DIR}/${FILIN} ${DATE1} ${DATE2} ${DT_BINS_H} ${MODE} $((DT_BINS_H*2/3))"
else
    CMD="${EXE2} ${RGPS_DIR}/${FILIN} ${DATE1} ${DATE2} ${DT_BINS_H} ${MODE}"
fi

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
