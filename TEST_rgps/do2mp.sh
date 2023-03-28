#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/generate_quad_mesh.py"

ijob=0
mkdir -p logs


csf="_${RESKM}km"
if [ ${RESKM} -eq 10 ]; then csf=""; fi



# Populating nc files we can use:
list_nc=`\ls nc/SELECTION_RGPS_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??${csf}.nc`
nbf=`echo ${list_nc} | wc -w`
echo " => ${nbf} files => ${nbf} batches!"


for ff in ${list_nc}; do

    fb=`basename ${ff}`
    echo
    echo " *** Doing file ${fb}"

    # Number of records inside netCDF file:
    Nr=`ncdump -h ${ff} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
    echo "   => it has ${Nr} time records"
    if [ ${Nr} -ne 2 ]; then
        echo "  ==> bad! PRESKMently we expect 2 records!"
        exit
    fi



    flog="`echo ${fb} | sed -e s/'.nc'/''/g`_${RESKM}km"

    ijob=$((ijob+1))

    CMD="${EXE} ${ff} 0,1 ${RESKM}"
    echo "    ==> will launch:"; echo "     ${CMD}"; echo
    ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
    #sleep 1
    echo

    if [ $((ijob%NJPAR)) -eq 0 ]; then
        echo "Waiting! (ijob = ${ijob})...."
        wait
        echo; echo
    fi

done

wait
