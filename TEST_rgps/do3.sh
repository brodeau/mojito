#!/bin/bash

# Deformation only

. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

mkdir -p logs


# Populating the batches available:
listQ=`\ls npz/Q-mesh_RGPS_S???_dt${DT_BINS_H}_${YEAR}????t0_${YEAR}????_*km.npz`

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
    list=`\ls npz/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_${YEAR}????t0_${YEAR}????_*km.npz`
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
        lst=( `\ls npz/Q-mesh_RGPS_${cbtch}_dt${DT_BINS_H}_${dr}_${YEAR}????_${RESKM}km.npz` )
        nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "

        if [ ${nf} -eq 2 ]; then

            fQ1=${lst[0]}
            fQ2=${lst[1]}
            echo " ==> will use:"
            echo " * ${fQ1}"
            echo " * ${fQ2}"

            cdt1=`echo ${fQ1} | cut -d_ -f5`
            cdt2=`echo ${fQ2} | cut -d_ -f5`
            flog="S${cbtch}_dt${DT_BINS_H}_${cdt1}-${cdt2}_${RESKM}km"

            if [ ${RESKM} -ge 50 ] && [ ${DT_BINS_H} -ge 72 ]; then
                #CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600*2/3))"
                CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600))"
            else
                CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600/2))"
            fi

            echo "  ==> ${CMD}"; echo
            ${CMD}
            echo; echo

        fi
        exit
    done

done
