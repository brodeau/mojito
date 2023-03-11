#!/bin/bash

# Deformation only

. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

LIST_RES="10" ; #fixme !!!

mkdir -p logs

ijob=0

# Populating the batches available:
listQ=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_S???_${YEAR}????t0_${YEAR}????_*km.npz`

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

for cres in ${LIST_RES}; do

    echo; echo
    echo " *** ${cres} km ***"
    echo

    for cbtch in ${list_btch}; do

        #  Q-mesh_RGPS_S000_19970104t0_19970104.npz
        list=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_${YEAR}????t0_${YEAR}????_${cres}km.npz`
        nbf=`echo ${list} | wc -w`

        echo " *** Number of files for Batch ${cbtch} = ${nbf}"

        list_date_ref=""
        for ff in ${list}; do
            date_ref=`echo ${ff} | cut -d_ -f6`
            list_date_ref+=" ${date_ref}"
        done
        list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

        echo; echo " *** List of reference dates for Batch${cbtch}:"; echo "${list_date_ref}"; echo

        for dr in ${list_date_ref}; do
            echo
            lst=( `\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${cbtch}_${dr}_${YEAR}????_${cres}km.npz` )
            nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
            
            if [ ${nf} -eq 2 ]; then

                fQ1=${lst[0]}
                fQ2=${lst[1]}
                echo " ==> will use:"
                echo " * ${fQ1}"
                echo " * ${fQ2}"
                
                flog=`basename ${fQ1}`
                flog=`echo ${flog} | sed -e s/".npz"/""/g`

                ijob=$((ijob+1))

                CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600/2))"
                echo "  ==> ${CMD}"; echo
                ${CMD} 1>logs/out_${flog}.out 2>logs/err_${flog}.err &
                echo; echo

                if [ $((ijob%NJPAR)) -eq 0 ]; then
                    echo "Waiting! (ijob = ${ijob})...."
                    wait
                    echo; echo
                fi

            fi

        done

    done

done

wait
