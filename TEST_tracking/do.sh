#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/ice_part_tracker.py"


CMD="${EXE} ${FSI3IN} ${FNMM}"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo

${CMD}
