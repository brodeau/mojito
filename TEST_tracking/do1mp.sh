#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/tracking/si3_part_tracker.py"

# 1/ populate the proper NC files to seed from:
dirin="${MOJITO_DIR}/TEST_brgps/nc"; echo ${dirin}

list_nc=`\ls ${dirin}/SELECTION_RGPS_*_${YEAR}????h??_${YEAR}????h??.nc`

nbf=`echo ${list_nc} | wc -w`

echo " *** We have ${nbf} files !"

mkdir -p ./logs

ijob=0

for fnc in ${list_nc}; do

    fb=`basename ${fnc}`
    echo "   * File: ${fb} :"

    # Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
    CMD="${EXE} ${FSI3IN} ${FNMM} ${fnc}" ; # with nc file for init seed...
    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo
    # SELECTION_RGPS_S000_19970104h15_19970107h15
    clog=`basename ${fnc} | sed -e s/'SELECTION_RGPS_'/''/g -e s/'.nc'/''/g`

    ${CMD} 1>./logs/out_${clog}.out 2>./logs/err_${clog}.err &

    ijob=$((ijob+1))

    sleep 3

    if [ $((ijob%NJPAR)) -eq 0 ]; then
        echo "Waiting! (ijob = ${ijob})...."
        wait
        echo; echo
    fi
    
done

wait
