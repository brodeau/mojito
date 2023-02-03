#!/bin/bash

# Deformation only

. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

mkdir -p logs

ijob=0

echo; echo; echo

cstr=`printf "%03d" ${istr}`

root_fn="Q-mesh_${NEMOCONF}_${EXPRMNT}_HSS${HSS}_TSS${TSS}"

# Q-mesh_NANUK4_ICE-BBM00_1h_HSS5_TSS72_19970107t0_19970107.npz

#list=( `\ls npz/Q-mesh_S${cstr}*_*_${dd}km_*.npz` )

list=( `ls npz/${root_fn}_${YEAR}????t0_${YEAR}????.npz` )

nbf=`echo ${list[*]} | wc -w`
echo " *** Number of files: ${nbf}"

nbf=`expr ${nbf} - 1`


if [ "${nbf}" != "" ]; then
    list_date_t0=""
    for ff in ${list[*]}; do
        date_t0=`echo ${ff} | cut -d_ -f7` ; #        echo "date_t0=${date_t0}"
        list_date_t0+=" ${date_t0}"
    done
    list_date_t0=$(echo ${list_date_t0} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

    echo; echo " *** List of reference dates for Stream${cstr}:"; echo "${list_date_t0}"; echo
    
    for dr in ${list_date_t0}; do
        echo
        lst=(`\ls npz/${root_fn}_${dr}_${YEAR}????.npz`) ; # echo ${lst[*]}
        nf=`echo ${lst[*]} | wc -w` ; echo " => ${nf} files "
        
        if [ ${nf} -eq 2 ]; then

            fQ1=${lst[0]}
            fQ2=${lst[1]}
            echo " ==> will use:"; echo "     * ${fQ1}"; echo "     * ${fQ2}"

            flog=`basename ${fQ1}`
            flog=`echo ${flog} | sed -e s/".npz"/""/g`

            ijob=$((ijob+1))

            CMD="${EXE} ${fQ1} ${fQ2} 0 ${MARKER_SIZE}"
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

fi
