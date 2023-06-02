#!/bin/bash

. ./conf.bash

EXE="python3 -u ${MOJITO_DIR}/generate_quad_mesh.py"
EXE2="python3 -u ${MOJITO_DIR}/generate_quad_rccl.py"

NSS=72 ; # Because we save every hourly time-steps in the netCDF files

cxtraRES=""
if [ "${LIST_RD_SS}" = "" ]; then
    cxtraRES="_${RESKM}km"
else
    cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
    cxtraRES="_${cr1}-${RESKM}km"
fi

mkdir -p ./logs

ijob=0

for NEMO_EXP in ${LIST_NEMO_EXP}; do

    # Populating nc files we can use:
    list_nc=`\ls nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking12_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??${cxtraRES}.nc`
    nbf=`echo ${list_nc} | wc -w`
    echo " => ${nbf} files => ${nbf} batches!"

    if [ "${ISEED_BASE}" = "quads" ] || [ "${ISEED_BASE}" = "defs" ]; then
        lRCCL=true
        EXE=${EXE2}
        # Populating npz files to get Quads identification from:
        if [ "${ISEED_BASE}" = "quads" ]; then
            list_seed_qnpz=`\ls ${DIRIN_PREPARED_RGPS}/npz/Q-mesh_RGPS_S???_dt${DT_BINS_H}_${YEAR}????-??h??t0_${YEAR}????-??h??${cxtraRES}${XTRASFX}.npz`
            # need to keep only the first occurence of each:
            lnew="" ; icpt=0
            for ff in ${list_seed_qnpz}; do
                if [ $((icpt%2)) -eq 0 ]; then lnew+="${ff} "; fi
                icpt=$((icpt+1))
            done
            list_seed_qnpz=${lnew}

        elif [ "${ISEED_BASE}" = "defs" ]; then
            list_seed_qnpz=`\ls ${DIRIN_PREPARED_RGPS}/npz/QUADSofDEF_RGPS_S???_dt${DT_BINS_H}_${YEAR}????${cxtraRES}${XTRASFX}.npz`
        fi
        #
        nbq=`echo ${list_seed_qnpz} | wc -w`
        echo " *** We have ${nbq} npz files to identify quadrangles:"
        echo ${list_seed_qnpz}
        echo
        if [ ${nbq} -ne ${nbf} ]; then
            echo "ERROR: that is not the same number as nc files!!!"; exit
        fi

        list_seed_qnpz=( ${list_seed_qnpz} )
    fi


    icpt=0
    for ff in ${list_nc}; do

        if ${lRCCL}; then fnpz=${list_seed_qnpz[${icpt}]}; fi

        if [ "${LIST_RD_SS}" = "" ]; then
            if [ -f ${ff} ]; then
                fb=`basename ${ff}`; echo; echo " *** Doing file ${fb}"
                ## Number of records inside netCDF file:
                #Nr=`ncdump -h ${ff} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
                #echo "   => it has ${Nr} time records"
                #if [ ${Nr} -le 65 ] || [ $((Nr/NSS)) -gt 1 ]; then
                #    echo "  ==> bad! Presently we expect ABOUT $((NSS+1)) records!"; exit
                #fi
                #lstrec="0,$((Nr-1))"
                lstrec="0,1"
                flog="quadgener_`echo ${fb} | sed -e s/'.nc'/''/g | sed -e s/"NEMO-SI3_${NEMO_CONF}_"/""/g`_${RESKM}km"
                ijob=$((ijob+1))
                if ${lRCCL}; then
                    CMD="${EXE} ${ff} ${lstrec} ${fnpz} ${RESKM} ${MODE}"
                else
                    CMD="${EXE} ${ff} ${lstrec} ${RESKM} ${MODE}"
                fi
                echo "    ==> will launch:"; echo "     ${CMD}"; echo
                ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
                echo
                if [ $((ijob%NJPAR)) -eq 0 ]; then
                    echo "Waiting! (ijob = ${ijob})...."
                    wait; echo; echo
                fi
            fi
        else
            for rdss in ${LIST_RD_SS}; do
                #fn=`echo ${ff} | sed -e "s|_${cr1}-${RESKM}km|_${rdss}-${RESKM}km|g"` ; # hour can be different!
                fn=`echo ${ff} | cut -d'_' -f1-7`
                fn=`basename ${fn}`
                fn=`\ls ./nc/${fn}_*_${rdss}-${RESKM}km.nc`
                if [ `echo ${fn} | wc -w` -ne 1 ]; then
                    echo "ERROR: less or more than 1 candidate for fn !!!"; exit
                fi                
                fnpzR=`echo ${fnpz} | sed -e "s|_${cr1}-${RESKM}km|_${rdss}-${RESKM}km|g"`
                #
                if [ -f ${fn} ]; then
                    fb=`basename ${fn}`; echo; echo " *** Doing file ${fb}"
                    lstrec="0,1"
                    flog="quadgener_`echo ${fb} | sed -e s/'.nc'/''/g | sed -e s/"NEMO-SI3_${NEMO_CONF}_"/""/g`_${RESKM}km"
                    ijob=$((ijob+1))
                    if ${lRCCL}; then
                        CMD="${EXE} ${fn} ${lstrec} ${fnpzR} ${RESKM} ${MODE} ${rdss}"
                    else
                        CMD="${EXE} ${fn} ${lstrec} ${RESKM} ${MODE} ${rdss}"
                    fi
                    echo "    ==> will launch:"; echo "     ${CMD}"; echo
                    ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
                    echo
                    if [ $((ijob%NJPAR)) -eq 0 ]; then
                        echo "Waiting! (ijob = ${ijob})...."
                        wait; echo; echo
                    fi
                else
                    echo "ERROR: file ${fn} is missing!!!"; exit
                fi
            done
        fi

        icpt=`expr ${icpt} + 1`

    done

done
wait
