#!/bin/bash

. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="python3 -u ${MOJITO_DIR}/generate_quad_mesh.py"

NSS=72 ; # Because we save every hourly time-steps in the netCDF files

cxtraRES=""
if [ ${RESKM} -gt 10 ]; then cxtraRES="_${RESKM}km"; fi

mkdir -p ./logs

ijob=0

for NEMO_EXP in ${LIST_NEMO_EXP}; do

    # Populating nc files we can use:
    list_nc=`\ls nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??${cxtraRES}.nc`
    nbf=`echo ${list_nc} | wc -w`
    echo " => ${nbf} files => ${nbf} batches!"


    for ff in ${list_nc}; do

        fb=`basename ${ff}`
        echo
        echo " *** Doing file ${fb}"

        # Number of records inside netCDF file:
        Nr=`ncdump -h ${ff} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
        echo "   => it has ${Nr} time records"
        if [ ${Nr} -le 65 ] || [ $((Nr/NSS)) -gt 1 ]; then
            echo "  ==> bad! Presently we expect ABOUT $((NSS+1)) records!"
            exit
        fi

        lstrec="0,$((Nr-1))"

        flog="quadgener_`echo ${fb} | sed -e s/'.nc'/''/g | sed -e s/"NEMO-SI3_${NEMO_CONF}_"/""/g`_${RESKM}km"

        ijob=$((ijob+1))

        CMD="${EXE} ${ff} ${lstrec} ${RESKM}"
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

done
wait
