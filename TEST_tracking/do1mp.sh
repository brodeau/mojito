#!/bin/sh

. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="${MOJITO_DIR}/tracking/si3_part_tracker.py"

# 1/ populate the proper NC files to seed from:
echo " * Will get RGPS seeding info in: ${DIRIN_PREPARED_RGPS} for RESKM = ${RESKM}"
cxtraRES=""
if [ ${RESKM} -gt 10 ]; then cxtraRES="_${RESKM}km"; fi

list_nc=`\ls ${DIRIN_PREPARED_RGPS}/SELECTION_RGPS_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??${cxtraRES}.nc`

nbf=`echo ${list_nc} | wc -w`

echo " *** We have ${nbf} files !"

echo ${list_nc}

mkdir -p ./logs

ijob=0

for fnc in ${list_nc}; do

    fb=`basename ${fnc}`
    echo "   * File: ${fb} :"

    # Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
    CMD="${EXE} ${FSI3IN} ${FNMM} ${fnc}" ; # with nc file for init seed...
    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo

    clog=`basename ${fnc} | sed -e s/"SELECTION_RGPS_"/"${NEMO_EXP}_"/g -e s/".nc"/""/g`

    ${CMD} 1>./logs/out_${clog}.out 2>./logs/err_${clog}.err &

    ijob=$((ijob+1))

    sleep 1

    if [ $((ijob%NJPAR)) -eq 0 ]; then
        echo "Waiting! (ijob = ${ijob})...."
        wait
        echo; echo
    fi
    
done

wait
