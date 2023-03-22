#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/statistics.py"

CMD="${EXE} npz ${DT_BINS_H} ${RESKM}"
echo "  ==> ${CMD}"; echo
${CMD}
echo; echo
