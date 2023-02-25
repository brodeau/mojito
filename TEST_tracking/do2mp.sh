#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/generate_quad_mesh.py"

LIST_RES=10; #fixme!

NSS=72 ; # Because we save every hourly time-steps in the netCDF files

# Populating nc files we can use:
list_nc=`\ls nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking_S???_${YEAR}????h??_${YEAR}????h??.nc`
nbf=`echo ${list_nc} | wc -w`
echo " => ${nbf} files => ${nbf} streams!"

ijob=0
mkdir -p logs

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

    for res in ${LIST_RES}; do

        flog="`echo ${fb} | sed -e s/'.nc'/''/g`_${res}km"

        ijob=$((ijob+1))

        CMD="${EXE} ${ff} ${lstrec} ${res}"
        echo "    ==> will launch:"; echo "     ${CMD}"; echo
        ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
        sleep 5
        echo
        
        if [ $((ijob%NJPAR)) -eq 0 ]; then
            echo "Waiting! (ijob = ${ijob})...."
            wait
            echo; echo
        fi
        
    done
    
done

wait





