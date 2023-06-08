#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/gatherDef.py"

CMD="${EXE} ./npz ${DT_BINS_H} ${RESKM} RGPS S"
echo "  ==> ${CMD}"; echo
${CMD}
echo; echo
