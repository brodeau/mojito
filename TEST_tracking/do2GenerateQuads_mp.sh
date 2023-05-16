#!/bin/bash

. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="python3 -u ${MOJITO_DIR}/generate_quad_mesh.py"

NSS=72 ; # Because we save every hourly time-steps in the netCDF files

cxtraRES=""
if [ "${LIST_RD_SS}" = "" ]; then
    if [ ${RESKM} -gt 10 ]; then cxtraRES="_${RESKM}km"; fi
else
    cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
    cxtraRES="_${cr1}-${RESKM}km"
fi

mkdir -p ./logs

ijob=0

for NEMO_EXP in ${LIST_NEMO_EXP}; do

    # Populating nc files we can use:
    list_nc=`\ls nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??${cxtraRES}.nc`
    nbf=`echo ${list_nc} | wc -w`
    echo " => ${nbf} files => ${nbf} batches!"


    for ff in ${list_nc}; do

        if [ "${LIST_RD_SS}" = "" ]; then
            if [ -f ${ff} ]; then
                fb=`basename ${ff}`; echo; echo " *** Doing file ${fb}"
                # Number of records inside netCDF file:
                Nr=`ncdump -h ${ff} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
                #if [ "${Nr}" != "" ]; then
                echo "   => it has ${Nr} time records"
                if [ ${Nr} -le 65 ] || [ $((Nr/NSS)) -gt 1 ]; then
                    echo "  ==> bad! Presently we expect ABOUT $((NSS+1)) records!"; exit
                fi
                lstrec="0,$((Nr-1))"
                flog="quadgener_`echo ${fb} | sed -e s/'.nc'/''/g | sed -e s/"NEMO-SI3_${NEMO_CONF}_"/""/g`_${RESKM}km"
                ijob=$((ijob+1))
                CMD="${EXE} ${ff} ${lstrec} ${RESKM} ${MODE}"
                echo "    ==> will launch:"; echo "     ${CMD}"; echo
                ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
                echo
                if [ $((ijob%NJPAR)) -eq 0 ]; then
                    echo "Waiting! (ijob = ${ijob})...."
                    wait; echo; echo
                fi
            fi
        else
            for rdss in ${LIST_RD_SS}; do
                fn=`echo ${ff} | sed -e "s|_${cr1}-${RESKM}km|_${rdss}-${RESKM}km|g"`
                if [ -f ${fn} ]; then
                    fb=`basename ${fn}`; echo; echo " *** Doing file ${fb}"
                    # Number of records inside netCDF file:
                    Nr=`ncdump -h ${fn} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
                    #if [ "${Nr}" != "" ]; then
                    echo "   => it has ${Nr} time records"
                    if [ ${Nr} -le 65 ] || [ $((Nr/NSS)) -gt 1 ]; then
                        echo "  ==> bad! Presently we expect ABOUT $((NSS+1)) records!"; exit
                    fi
                    lstrec="0,$((Nr-1))"
                    flog="quadgener_`echo ${fb} | sed -e s/'.nc'/''/g | sed -e s/"NEMO-SI3_${NEMO_CONF}_"/""/g`_${RESKM}km"
                    ijob=$((ijob+1))
                    CMD="${EXE} ${fn} ${lstrec} ${RESKM} ${MODE} ${rdss}"
                    echo "    ==> will launch:"; echo "     ${CMD}"; echo
                    ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
                    echo
                    if [ $((ijob%NJPAR)) -eq 0 ]; then
                        echo "Waiting! (ijob = ${ijob})...."
                        wait; echo; echo
                    fi
                fi
            done
        fi



    done

done
wait
