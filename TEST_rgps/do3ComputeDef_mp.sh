#!/bin/bash

# Deformation only

. ./conf.bash

EXE="python3 -u ${MOJITO_DIR}/deformation.py"
EXE2="python3 -u ${MOJITO_DIR}/gatherSlctBatch.py"

mkdir -p ./logs

ijob=0


csf="_${RESKM}km"
if [ "${LIST_RD_SS}" != "" ]; then
    cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
    csf="_${cr1}-${RESKM}km"
fi

# Populating the batches available:
listQ=`\ls msh/Q-mesh_RGPS_S???_dt${DT_BINS_H}_199[67]????t0_199[67]????${csf}.npz`

echo "${listQ}"

list_btch=""
for ff in ${listQ}; do
    list_btch+="`echo ${ff} | cut -d'_' -f3` "
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
        list=`\ls msh/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_199[67]????t0_199[67]????${csf}.npz 2>/dev/null`
        if [ "${list}" != "" ]; then

            nbf=`echo ${list} | wc -w`
            echo " *** Number of files for Batch ${cbtch} = ${nbf}"

            list_date_ref=""
            for ff in ${list}; do
                date_ref=`echo ${ff} | cut -d_ -f5`
                list_date_ref+=" ${date_ref}"
            done
            #list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !
            list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -u) ; # unique !

            echo; echo " *** List of reference dates for Batch${cbtch}:"; echo "${list_date_ref}"; echo


            for dr in ${list_date_ref}; do
                echo

                if [ "${LIST_RD_SS}" = "" ]; then
                    lst=( `\ls msh/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_${dr}_199[67]????${csf}.npz 2>/dev/null` )
                    if [ "${list}" != "" ]; then
                        nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
                        #
                        if [ ${nf} -eq 2 ]; then
                            fQ1=${lst[0]}; fQ2=${lst[1]}
                            echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                            cdt1=`echo ${fQ1} | cut -d_ -f5`; cdt2=`echo ${fQ2} | cut -d_ -f5`

                            flog="def_S${cbtch}_dt${DT_BINS_H}_${cdt1}-${cdt2}_${RESKM}km"
                            #
                            if [ ${RESKM} -ge 50 ] && [ ${DT_BINS_H} -ge 72 ]; then
                                CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600))         ${MODE} ${DEF_EXPORT}"
                            else
                                CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600/2*20/15)) ${MODE} ${DEF_EXPORT}"
                            fi
                            echo "  ==> ${CMD}"; echo
                            ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
                            echo; echo
                            ijob=$((ijob+1))
                            if [ $((ijob%NJPAR)) -eq 0 ]; then
                                echo "Waiting! (ijob = ${ijob})...."
                                wait; echo; echo
                            fi
                        fi ; # if [ ${nf} -eq 2 ]
                    fi
                    #
                else
                    #
                    for rdss in ${LIST_RD_SS}; do
                        lst=( `\ls msh/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_${dr}_199[67]????_${rdss}-${RESKM}km.npz 2>/dev/null` )
                        if [ "${list}" != "" ]; then
                            nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
                            #
                            if [ ${nf} -eq 2 ]; then
                                fQ1=${lst[0]}; fQ2=${lst[1]}
                                echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                                cdt1=`echo ${fQ1} | cut -d_ -f5`; cdt2=`echo ${fQ2} | cut -d_ -f5`

                                flog="def_S${cbtch}_dt${DT_BINS_H}_${cdt1}-${cdt2}_${rdss}-${RESKM}km"
                                ijob=$((ijob+1))
                                #
                                if [ ${RESKM} -ge 50 ] && [ ${DT_BINS_H} -ge 72 ]; then
                                    CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600))         ${MODE} ${DEF_EXPORT}"
                                else
                                    CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600/2*20/15)) ${MODE} ${DEF_EXPORT}"
                                fi
                                echo "  ==> ${CMD}"; echo
                                ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
                                echo; echo
                                if [ $((ijob%NJPAR)) -eq 0 ]; then
                                    echo "Waiting! (ijob = ${ijob})...."
                                    wait; echo; echo
                                fi
                            fi ; # if [ ${nf} -eq 2 ]
                        fi
                        #
                    done
                    #
                fi # if [ "${LIST_RD_SS}" = "" ]

            done

        fi

    done

fi ; # if [ "$1" != "2" ]


wait

if [ "${LIST_RD_SS}" != "" ] && [ "${USE_S}" != "" ]; then

    for cbtch in ${list_btch}; do



        CMD="${EXE2} ./npz ${cbtch} dt${DT_BINS_H} ${RESKM} RGPS"

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


fi

wait
