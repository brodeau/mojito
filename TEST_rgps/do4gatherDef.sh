#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/gatherDef.py"

CMD="${EXE} ./npz ${DT_BINS_H} ${RESKM} RGPS"
echo "  ==> ${CMD}"; echo
${CMD}
echo; echo
