#!/bin/bash

# Deformation only

. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

mkdir -p logs

ijob=0

for dd in ${LIST_RES}; do

    echo; echo
    echo " *** ${dd} km ***"
    echo

    for st in ${LIST_STREAM}; do
        echo; echo; echo

        list=( `\ls npz/Q-mesh_S${st}*_*_${dd}km_*.npz` )
        nbf=`echo ${list[*]} | wc -w`
        nbf=`expr ${nbf} - 1`

        echo " *** Number of files for Stream ${st} = ${nbf}"

        if [ "${nbf}" != "" ]; then
            list_date_ref=""            
            for ff in ${list[*]}; do
                date_ref=`echo ${ff} | cut -d_ -f3`
                list_date_ref+=" ${date_ref}"
            done
            list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

            echo; echo " *** List of reference dates for Stream${st}:"; echo "${list_date_ref}"; echo

            for dr in ${list_date_ref}; do
                echo
                lst=(`\ls npz/Q-mesh_S${st}_${dr}_*_${dd}km_*.npz`) ; # echo ${lst[*]}
                nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
                if [ ${nf} -eq 2 ]; then

                    fQ1=${lst[0]}
                    fQ2=${lst[1]}
                    echo " ==> will use:"; echo "     * ${fQ1}"; echo "     * ${fQ2}"

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

        fi

    done

done
