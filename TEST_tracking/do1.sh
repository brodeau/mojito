#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/ice_part_tracker.py"

# 1/ populate the proper NC files to seed from:
dirin="${MOJITO_DIR}/TEST_brgps/nc"; echo ${dirin}

list_nc=`\ls ${dirin}/SELECTION_buoys_RGPS_*_${YEAR}????h??_${YEAR}????h??.nc`

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
