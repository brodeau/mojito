#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/ice_part_tracker.py"


if [ "$1" = "d" ]; then
    CMD="${EXE} ${FSI3IN} ${FNMM}"
else
    CMD="${EXE} ${FSI3IN} ${FNMM} ${FNCSEED}" ; # with nc file for init seed...
fi

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo

${CMD}
