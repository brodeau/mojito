#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/statistics.py"

CMD="${EXE} npz RGPS"
echo "  ==> ${CMD}"; echo
${CMD}
echo; echo
