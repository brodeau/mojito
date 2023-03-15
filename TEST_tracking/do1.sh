#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/tracking/si3_part_tracker.py"

# 1/ populate the proper NC files to seed from:
echo " * Will get RGPS seeding info in: ${DIRIN_PREPARED_RGPS}"

list_nc=`\ls ${DIRIN_PREPARED_RGPS}/SELECTION_RGPS_*_${YEAR}????h??_${YEAR}????h??.nc`

nbf=`echo ${list_nc} | wc -w`

echo " *** We have ${nbf} files !"

for fnc in ${list_nc}; do

    fb=`basename ${fnc}`
    echo "   * File: ${fb} :"

    # Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
    CMD="${EXE} ${FSI3IN} ${FNMM} ${fnc}" ; # with nc file for init seed...
    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo

    ${CMD}

done

wait
