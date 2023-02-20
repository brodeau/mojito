#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/statistics.py"

#if [ "$1" = "" ]; then
#    echo "USAGE: $0 <file_pos_mojito.nc>"
#    exit
#fi
#FIN="$1"

CMD="${EXE} ./npz ${NEMO_EXP}"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo

${CMD}
