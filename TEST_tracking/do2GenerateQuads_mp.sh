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

# dt${DT_BINS_H} => 1h

mkdir -p ./logs

ijob=0

for NEMO_EXP in ${LIST_NEMO_EXP}; do

    # Populating nc files we can use:
    list_nc=`\ls nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking12_S???_dt${DT_BINS_H}_1h_199[67]????h??_199[67]????h??${cxtraRES}.nc`
    nbf=`echo ${list_nc} | wc -w`
    echo " => ${nbf} nc files => ${nbf} batches!"

    if [ "${ISEED_BASE}" = "quads" ] || [ "${ISEED_BASE}" = "defs" ]; then
        lRCCL=true
        EXE=${EXE2}
        # Populating npz files to get Quads identification from:
        if [ "${ISEED_BASE}" = "quads" ]; then
            list_seed_qnpz=`\ls ${DIRIN_PREPARED_RGPS}/msh/Q-mesh_RGPS_S???_dt${DT_BINS_H}_199[67]????-??h??t0_199[67]????-??h??${cxtraRES}${XTRASFX}.npz`
            # need to keep only the first occurence of each:
            lnew="" ; icpt=0
            for ff in ${list_seed_qnpz}; do
                if [ $((icpt%2)) -eq 0 ]; then lnew+="${ff} "; fi
                icpt=$((icpt+1))
            done
            list_seed_qnpz=${lnew}

        elif [ "${ISEED_BASE}" = "defs" ]; then
            echo "LOLO2!"
            list_seed_qnpz=`\ls ${DIRIN_PREPARED_RGPS}/npz/QUADSofDEF_RGPS_S???_dt${DT_BINS_H}_199[67]????${cxtraRES}${XTRASFX}.npz`
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
        echo; echo " *** Based on file ${ff} !"
        SB=`basename ${ff} | cut -d'_' -f 5`
        if [ "`echo ${SB} | cut -c1-1`" != "S" ]; then
            echo " ERROR: SB does not seem to get the batch string!"; exit
        fi
        echo "      ==> batch is: ${SB}"

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
                #
                if [ "${ISEED_BASE}" = "defs" ]; then
                    fnc=`\ls nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking12_${SB}_dt${DT_BINS_H}_1h_199[67]????h??_199[67]????h??_${rdss}-${RESKM}km.nc 2>/dev/null`
                    fnp=`\ls ${DIRIN_PREPARED_RGPS}/npz/QUADSofDEF_RGPS_${SB}_dt${DT_BINS_H}_199[67]????_${rdss}-${RESKM}km.npz 2>/dev/null`
                else
                    fnc=`echo ${ff} | sed -e "s|_${cr1}-${RESKM}km|_${rdss}-${RESKM}km|g"`
                    fnp=`echo ${fnpz} | sed -e "s|_${cr1}-${RESKM}km|_${rdss}-${RESKM}km|g"`
                fi
                #
                if [ "${fnc}" != "" ] && [ -f "${fnc}" ]; then
                    #
                    if [ "${fnp}" != "" ] && [ -f "${fnp}" ]; then
                        #
                        if [ `echo ${fnc} | wc -w` -ne 1 ]; then
                            echo "ERROR: less or more than 1 candidate for fnc !!!"; echo "  ==> ${fnc}"; exit
                        fi

                        #
                        fb=`basename ${fnc}`; echo; echo " *** Doing file ${fb}"
                        lstrec="0,1"
                        flog="quadgener_`echo ${fb} | sed -e s/'.nc'/''/g | sed -e s/"NEMO-SI3_${NEMO_CONF}_"/""/g`_${RESKM}km"
                        ijob=$((ijob+1))
                        if ${lRCCL}; then
                            CMD="${EXE} ${fnc} ${lstrec} ${fnp} ${RESKM} ${MODE} ${rdss}"
                        else
                            CMD="${EXE} ${fnc} ${lstrec} ${RESKM} ${MODE} ${rdss}"
                        fi
                        echo "    ==> will launch:"; echo "     ${CMD}"; echo ; #exit;#lolo
                        ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
                        echo
                        if [ $((ijob%NJPAR)) -eq 0 ]; then
                            echo "Waiting! (ijob = ${ijob})...."
                            wait; echo; echo
                        fi
                    else
                        echo "WARNING: npz file *_${SB}_*_${rdss}-${RESKM}km.npz was not there..."
                    fi
                else
                    echo "WARNING: nc file *_${SB}_*_${rdss}-${RESKM}km.nc was not there..."
                fi
            done
        fi

        icpt=`expr ${icpt} + 1`

    done

done
wait
