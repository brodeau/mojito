#!/bin/bash

# Deformation only

. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

mkdir -p ./logs

ijob=0


csf="_${RESKM}km"
if [ "${LIST_RD_SS}" != "" ]; then
    cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
    csf="_${cr1}-${RESKM}km"
fi

# Populating the batches available:
listQ=`\ls npz/Q-mesh_RGPS_S???_dt${DT_BINS_H}_${YEAR}????t0_${YEAR}????${csf}.npz`

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

for cbtch in ${list_btch}; do

    #  Q-mesh_RGPS_S000_19970104t0_19970104.npz
    list=`\ls npz/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_${YEAR}????t0_${YEAR}????${csf}.npz`
    nbf=`echo ${list} | wc -w`

    echo " *** Number of files for Batch ${cbtch} = ${nbf}"

    list_date_ref=""
    for ff in ${list}; do
        date_ref=`echo ${ff} | cut -d_ -f5`
        list_date_ref+=" ${date_ref}"
    done
    list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

    echo; echo " *** List of reference dates for Batch${cbtch}:"; echo "${list_date_ref}"; echo

    for dr in ${list_date_ref}; do
        echo

        if [ "${LIST_RD_SS}" = "" ]; then
            lst=( `\ls npz/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_${dr}_${YEAR}????${csf}.npz` )
            nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
            #
            if [ ${nf} -eq 2 ]; then
                fQ1=${lst[0]}; fQ2=${lst[1]}
                echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                cdt1=`echo ${fQ1} | cut -d_ -f5`; cdt2=`echo ${fQ2} | cut -d_ -f5`

                flog="def_S${cbtch}_dt${DT_BINS_H}_${cdt1}-${cdt2}_${RESKM}km"
                ijob=$((ijob+1))
                #
                if [ ${RESKM} -ge 50 ] && [ ${DT_BINS_H} -ge 72 ]; then
                    CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600))"
                else
                    CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600/2*20/19))"
                fi
                echo "  ==> ${CMD}"; echo
                ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
                echo; echo
                if [ $((ijob%NJPAR)) -eq 0 ]; then
                    echo "Waiting! (ijob = ${ijob})...."
                    wait; echo; echo
                fi
            fi ; # if [ ${nf} -eq 2 ]
            #
        else
            #
            for rdss in ${LIST_RD_SS}; do
                lst=( `\ls npz/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_${dr}_${YEAR}????_${rdss}-${RESKM}km.npz` )
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
                        CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600))"
                    else
                        CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600/2*20/19))"
                    fi
                    echo "  ==> ${CMD}"; echo
                    ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
                    echo; echo
                    if [ $((ijob%NJPAR)) -eq 0 ]; then
                        echo "Waiting! (ijob = ${ijob})...."
                        wait; echo; echo
                    fi
                fi ; # if [ ${nf} -eq 2 ]
                #
            done
            #
        fi # if [ "${LIST_RD_SS}" = "" ]

    done

done

wait
