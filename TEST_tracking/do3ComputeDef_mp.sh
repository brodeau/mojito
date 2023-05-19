#!/bin/bash

# Deformation only
. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

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
    listQ=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_S???_dt${DT_BINS_H}_${YEAR}????-??h??t0_${YEAR}????-??h??${csf}.npz`

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

    for cbtch in ${list_btch}; do

        #  Q-mesh_RGPS_S000_19970104t0_19970104.npz
        list=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_${YEAR}????-??h??t0_${YEAR}????-??h??${csf}.npz`
        nbf=`echo ${list} | wc -w`

        echo " *** Number of files for Batch ${cbtch} = ${nbf}"

        list_date_ref=""
        for ff in ${list}; do
            date_ref=`echo ${ff} | cut -d_ -f7`
            list_date_ref+=" ${date_ref}"
        done
        list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

        echo; echo " *** List of reference dates for Batch${cbtch}:"; echo "${list_date_ref}"; echo

        for dr in ${list_date_ref}; do
            echo

            if [ "${LIST_RD_SS}" == "" ]; then
                lst=( `\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_${dr}_${YEAR}????-??h??${csf}.npz 2>/dev/null` )
                if [ "${lst}" != "" ]; then
                    nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
                    #
                    if [ ${nf} -eq 2 ]; then
                        fQ1=${lst[0]}
                        fQ2=${lst[1]}
                        echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                        flog=`basename ${fQ1}`; flog=`echo ${flog} | sed -e s/".npz"/""/g`
                        ijob=$((ijob+1))
                        CMD="${EXE} ${fQ1} ${fQ2} 0 ${MODE}"
                        echo "  ==> ${CMD}"; echo
                        ${CMD} 1>logs/out_${flog}.out 2>logs/err_${flog}.err &
                        echo; echo
                        if [ $((ijob%NJPAR)) -eq 0 ]; then
                            echo "Waiting! (ijob = ${ijob})...."
                            wait; echo; echo
                        fi
                    fi
                fi
            else

                for rdss in ${LIST_RD_SS}; do
                    #
                    lst=( `\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_dt${DT_BINS_H}_${dr}_${YEAR}????-??h??_${rdss}-${RESKM}km.npz 2>/dev/null` )
                    if [ "${lst}" != "" ]; then
                        nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
                        #
                        if [ ${nf} -eq 2 ]; then
                            fQ1=${lst[0]}
                            fQ2=${lst[1]}
                            echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                            flog=`basename ${fQ1}`; flog=`echo ${flog} | sed -e s/".npz"/""/g`
                            ijob=$((ijob+1))
                            CMD="${EXE} ${fQ1} ${fQ2} 0"
                            echo "  ==> ${CMD}"; echo
                            ${CMD} 1>logs/out_${flog}.out 2>logs/err_${flog}.err &
                            echo; echo
                            if [ $((ijob%NJPAR)) -eq 0 ]; then
                                echo "Waiting! (ijob = ${ijob})...."
                                wait; echo; echo
                            fi
                        fi
                    fi

                done

            fi

        done

    done

done

wait
