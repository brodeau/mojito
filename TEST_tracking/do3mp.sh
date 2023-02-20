#!/bin/bash

# Deformation only

. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

mkdir -p logs

ijob=0


RESKM=10

cstr="NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}"

list=( `\ls npz/Q-mesh_${cstr}*.npz | grep ${YEAR}` )
nbf=`echo ${list[*]} | wc -w`
nbf=`expr ${nbf} - 1`


echo " *** Number of files = ${nbf}"


if [ "${nbf}" != "" ]; then
    list_date_ref=""
    for ff in ${list[*]}; do
        date_ref=`echo ${ff} | cut -d_ -f5`
        list_date_ref+=" ${date_ref}"
    done
    list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

    echo; echo " *** List of reference dates for ${cstr}:"; echo "${list_date_ref}"; echo

    
    for dr in ${list_date_ref}; do
        echo
        lst=(`\ls npz/Q-mesh_${cstr}_*${dr}_*.npz`) ;  #echo ${lst[*]}
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
