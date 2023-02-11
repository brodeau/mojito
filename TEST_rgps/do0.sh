#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/rgps/01_interp_to_calendar.py"

YEAR=`echo ${DATE1} | cut -c1-4`

#if [ `echo ${DATE2} | cut -c1-4` -ne ${YEAR} ]; then
#    echo "ERROR: DATE1 and DATE2 have 2 different years!" ; exit
#fi

MMDD1=`echo ${DATE1} | cut -c5-8`
MMDD2=`echo ${DATE2} | cut -c5-8`

CMD="${EXE} ${RGPS_DIR}/${FILIN} ${DATE1} ${DATE2} ${DT_BINS_H}"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
