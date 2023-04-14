#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/coarsify_point_cloud.py"

ijob=0
mkdir -p logs

# Populating nc files we can use:
list_nc=`\ls nc/SELECTION_RGPS_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??.nc`
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
        echo "  ==> bad! Presently we expect 2 records!"
        exit
    fi


    if [ "${LIST_RD_SS}" = "" ]; then
        fout=`echo ${ff} | sed -e "s|.nc|_${RESKM}km.nc|g"`
        if [ ! -f ${fout} ]; then
            flog="coarsify_`echo ${fb} | sed -e s/'.nc'/''/g`_${RESKM}km"
            CMD="${EXE} ${ff} ${RESKM}"
            ijob=$((ijob+1))
            echo "    ==> will launch:"; echo "     ${CMD}"; echo
            ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
            #
            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait
                echo; echo
            fi
        fi
    else
        for rdss in ${LIST_RD_SS}; do
            fout=`echo ${ff} | sed -e "s|.nc|_${rdss}-${RESKM}km.nc|g"`
            if [ ! -f ${fout} ]; then
                flog="coarsify_`echo ${fb} | sed -e s/'.nc'/''/g`_${rdss}-${RESKM}km"
                CMD="${EXE} ${ff} ${RESKM} ${rdss} "
                ijob=$((ijob+1))
                echo "    ==> will launch:"; echo "     ${CMD}"; echo
                ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
                #
                if [ $((ijob%NJPAR)) -eq 0 ]; then
                    echo "Waiting! (ijob = ${ijob})...."
                    wait
                    echo; echo
                fi
            fi
        done
    fi


done

wait
