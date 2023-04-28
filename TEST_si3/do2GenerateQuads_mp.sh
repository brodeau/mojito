#!/bin/bash

#. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="python3 -u ${MOJITO_DIR}/generate_quad_mesh.py"

NSS=72 ; # Because we save every hourly time-steps in the netCDF files




mkdir -p ./logs

ijob=0

for RESKM in ${LCOARSEN[*]}; do

    cxtraRES="_${RESKM}km"

    for NEMO_EXP in ${LIST_NEMO_EXP}; do

        # Populating nc files we can use:
        list_nc=`\ls nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking_nemoTsi3_idlSeed_${YEAR}????h??_${YEAR}????h??${cxtraRES}.nc`
        #list_nc=`\ls nc/NEMO-SI3_NANUK4_BBM00_tracking_nemoTsi3_idlSeed_19970107h00_19970110h00_640km.nc`
        nbf=`echo ${list_nc} | wc -w`
        echo " => ${nbf} files!"
        echo "${list_nc}"; echo


        for ff in ${list_nc}; do

            if [ -f ${ff} ]; then
                fb=`basename ${ff}`; echo; echo " *** Doing file ${fb}"
                # Number of records inside netCDF file:
                Nr=`ncdump -h ${ff} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
                #if [ "${Nr}" != "" ]; then
                echo "   => it has ${Nr} time records"
                if [ ${Nr} -le 65 ] || [ $((Nr/NSS)) -gt 1 ]; then
                    echo " WARNING:  ==> bad! Presently we expect ABOUT $((NSS+1)) records!";echo
                else
                    lstrec="0,$((Nr-1))"
                    flog="quadgener_`echo ${fb} | sed -e s/'.nc'/''/g | sed -e s/"NEMO-SI3_${NEMO_CONF}_"/""/g`_${RESKM}km"
                    ijob=$((ijob+1))
                    CMD="${EXE} ${ff} ${lstrec} ${RESKM}"
                    echo "    ==> will launch:"; echo "     ${CMD}"; echo

                    ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
                    echo
                    if [ $((ijob%NJPAR)) -eq 0 ]; then
                        echo "Waiting! (ijob = ${ijob})...."
                        wait; echo; echo
                    fi
                fi
            fi



        done

    done
done
wait
