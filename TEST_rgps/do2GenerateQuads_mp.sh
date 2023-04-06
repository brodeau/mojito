#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/generate_quad_mesh.py"

ijob=0
mkdir -p logs



csf="_${RESKM}km"
if [ "${LIST_RD_SS}" != "" ]; then
    cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
    csf="_${cr1}-${RESKM}km"
fi


if [ ${RESKM} -eq 10 ]; then csf=""; fi



# Populating nc files we can use:
list_nc=`\ls nc/SELECTION_RGPS_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??${csf}.nc`

nbf=`echo ${list_nc} | wc -w`
echo " => ${nbf} files => ${nbf} batches!"

echo; echo ${list_nc}; echo

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


    if [ "${LIST_RD_SS}" = "" ]; then
        flog="`echo ${fb} | sed -e s/'.nc'/''/g`_${RESKM}km"
        ijob=$((ijob+1))
        CMD="${EXE} ${ff} 0,1 ${RESKM}"
        echo "    ==> will launch:"; echo "     ${CMD}"; echo
        ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
        if [ $((ijob%NJPAR)) -eq 0 ]; then
            echo "Waiting! (ijob = ${ijob})...."
            wait; echo; echo
        fi
    else
        for rdss in ${LIST_RD_SS}; do
            flog="`echo ${fb} | sed -e s/'.nc'/''/g`_${rdss}-${RESKM}km"
            ijob=$((ijob+1))
            #
            fn=`echo ${ff} | sed -e "s|_${cr1}-${RESKM}km|_${rdss}-${RESKM}km|g"`            
            #
            CMD="${EXE} ${fn} 0,1 ${RESKM} ${rdss}"
            echo "    ==> will launch:"; echo "     ${CMD}"; echo
            ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait; echo; echo
            fi
        done
    fi

done

wait
