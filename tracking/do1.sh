#!/bin/sh

. ../TEST_tracking/conf.bash


DATE1=19970101
DATE2=19970131


EXE="${MOJITO_DIR}/tracking/si3_part_tracker.py"

# 1/ populate the proper NC files to seed from:
dirin="./nc"; echo ${dirin}

list_nc=`\ls ${dirin}/*.nc`

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
