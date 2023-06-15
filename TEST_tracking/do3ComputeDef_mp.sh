#!/bin/bash

# Deformation only

. ./conf.bash

EXE="python3 -u ${MOJITO_DIR}/deformation.py"
EXE2="python3 -u ${MOJITO_DIR}/gatherSlctBatch.py"


#MAX_T_DEV=7200 ; # 2 hours
MAX_T_DEV=$((DT_BINS_H*3600)) ; # like RGPS !!!

mkdir -p ./logs

ijob=0

for NEMO_EXP in ${LIST_NEMO_EXP}; do
    echo; echo

    csf="_${RESKM}km"
    if [ "${LIST_RD_SS}" != "" ]; then
        cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
        csf="_${cr1}-${RESKM}km"
    fi

    # Populating the batches available:
    listQ=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_S???_dt${DT_BINS_H}_199?????t0_199?????${csf}.npz`

    echo "${listQ}"

    list_btch=""
    for ff in ${listQ}; do
        list_btch+="`echo ${ff} | cut -d'_' -f5` "
    done

    # Removing double of occurences:
    list_btch=$(echo ${list_btch} | tr ' ' '\n' | sort -u)
    echo ${list_btch}
    nbs=`echo ${list_btch} | wc -w`
    echo " ==> ${nbs} batches!" ; echo

    echo; echo
    echo " *** ${RESKM} km ***"
    echo


    if [ "$1" != "2" ]; then

        for cbtch in ${list_btch}; do

            #  Q-mesh_RGPS_S000_19970104t0_19970104.npz
            list=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_199?????t0_199?????${csf}.npz`
            nbf=`echo ${list} | wc -w`

            echo " *** Number of files for Batch ${cbtch} with suffix ${csf}.npz = ${nbf}"

            #list_date_ref=""
            #for ff in ${list}; do
            #    date_ref=`echo ${ff} | cut -d_ -f7`
            #    list_date_ref+=" ${date_ref}"
            #done
            #list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

            #echo; echo " *** List of reference dates for Batch${cbtch}:"; echo "${list_date_ref}"; echo

            #exit



            #for dr in ${list_date_ref}; do
            #    echo

            if [ "${LIST_RD_SS}" == "" ]; then
                clst=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_199?????t0_199?????${csf}.npz 2>/dev/null`
                if [ "${clst}" != "" ]; then
                    lst=( ${clst} )
                    nbf=`echo ${lst[*]} | wc -w`
                    if [ ${nbf} -ne 2 ]; then echo "ERROR: we do not have 2 files!!!! => ${lst[*]}"; exit; fi
                    #
                    fQ1=${lst[0]}
                    fQ2=${lst[1]}
                    echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                    flog="def__`basename ${fQ1}`"; flog=`echo ${flog} | sed -e s/".npz"/""/g`
                    ijob=$((ijob+1))
                    CMD="${EXE} ${fQ1} ${fQ2} ${MAX_T_DEV} ${MODE}"
                    echo "  ==> ${CMD}"; echo
                    ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
                    echo; echo
                    if [ $((ijob%NJPAR)) -eq 0 ]; then
                        echo "Waiting! (ijob = ${ijob})...."
                        wait; echo; echo
                    fi
                else
                    echo "WARNING: No files found for: Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_199?????t0_199?????_${rdss}-${RESKM}km.npz !"
                    echo
                fi
            else

                for rdss in ${LIST_RD_SS}; do
                    #
                    clst=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_199?????t0_199?????_${rdss}-${RESKM}km.npz 2>/dev/null`
                    if [ "${clst}" != "" ]; then
                        lst=( ${clst} )
                        nbf=`echo ${lst[*]} | wc -w`
                        if [ ${nbf} -ne 2 ]; then echo "ERROR: we do not have 2 files!!!! => ${lst[*]}"; exit; fi

                        if [ "${lst}" != "" ]; then
                            fQ1=${lst[0]}
                            fQ2=${lst[1]}
                            echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                            flog="def__`basename ${fQ1}`"; flog=`echo ${flog} | sed -e s/".npz"/""/g`
                            ijob=$((ijob+1))
                            CMD="${EXE} ${fQ1} ${fQ2} ${MAX_T_DEV} ${MODE}"
                            echo "  ==> ${CMD}"; echo ; #exit;#lolo
                            ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
                            echo; echo
                            if [ $((ijob%NJPAR)) -eq 0 ]; then
                                echo "Waiting! (ijob = ${ijob})...."
                                wait; echo; echo
                            fi
                        fi
                    else
                        echo "WARNING: No files found for: Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_199?????t0_199?????_${rdss}-${RESKM}km.npz !"
                        echo
                    fi

                done

            fi

            #done

        done

    fi ; # if [ "$1" != "2" ]


done

wait




if [ "${LIST_RD_SS}" != "" ] && [ "${USE_S}" != "" ]; then

    for NEMO_EXP in ${LIST_NEMO_EXP}; do
        echo; echo



        for cbtch in ${list_btch}; do



            CMD="${EXE2} ./npz ${cbtch} dt${DT_BINS_H} ${RESKM} NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}"

            flog="gatherSlctBatch_${cbtch}_dt${DT_BINS_H}_${RESKM}"

            echo "  ==> ${CMD}"; echo
            ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
            echo; echo
            ijob=$((ijob+1))
            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait; echo; echo
            fi

        done





    done
    wait

fi
