#!/bin/bash

. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="${MOJITO_DIR}/statistics.py"

#if [ "$1" = "" ]; then
#    echo "USAGE: $0 <file_pos_mojito.nc>"
#    exit
#fi
#FIN="$1"

CMD="${EXE} ./npz ${DT_BINS_H} ${RESKM} NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo

${CMD}
