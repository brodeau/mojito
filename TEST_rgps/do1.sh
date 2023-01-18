#!/bin/bash

. ./conf.bash

idebug=0

EXE="${MOJITO_DIR}/rgps/01_selection_xy.py"
EXEdbg="${MOJITO_DIR}/rgps/d01_follow_set_of_buoys.py"

YEAR=`echo ${DATE1} | cut -c1-4`

if [ `echo ${DATE2} | cut -c1-4` -ne ${YEAR} ]; then
    echo "ERROR: DATE1 and DATE2 have 2 different years!" ; exit
fi

MMDD1=`echo ${DATE1} | cut -c5-8`
MMDD2=`echo ${DATE2} | cut -c5-8`


if [ ${idebug} -eq 1 ]; then
    CMDdbg="${EXEdbg} ${RGPS_DIR}/${FILIN} ${YEAR} ${MMDD1} ${MMDD2}"
    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo
    ${CMDdbg}
    echo
    exit
fi



CMD="${EXE} ${RGPS_DIR}/${FILIN} ${YEAR} ${MMDD1} ${MMDD2} ${DT_BINS_H}"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
