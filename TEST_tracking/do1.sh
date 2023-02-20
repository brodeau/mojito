#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/ice_part_tracker.py"


if [ ! -f ${FNCSEED} ]; then
    echo " ERROR: file ${FNCSEED} is missing !"; exit
fi


if [ "$1" = "d" ]; then
    CMD="${EXE} ${FSI3IN} ${FNMM} ${DATE2}"
else
    CMD="${EXE} ${FSI3IN} ${FNMM} ${DATE2} ${FNCSEED}" ; # with nc file for init seed...
fi

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo

${CMD}
