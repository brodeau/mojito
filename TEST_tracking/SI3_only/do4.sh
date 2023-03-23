#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/statistics.py"


PRFX="DEFORMATIONS_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemoTsi3_${YEAR}"

#if [ "$1" = "" ]; then
#    echo "USAGE: $0 <file_pos_mojito.nc>"
#    exit
#fi
#FIN="$1"
#${NEMO_EXP}

CMD="${EXE} ./npz ${PRFX}"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
