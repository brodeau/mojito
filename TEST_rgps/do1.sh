#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/rgps/01_selection_xy.py"

CMD="${EXE} ${RGPS_DIR}/${FILIN} ${YEAR}"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
