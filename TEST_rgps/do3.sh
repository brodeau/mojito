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

        list=( `\ls npz/Q-mesh_S${st}*_*_${dd}km_*.npz` )
        nbf=`echo ${list[*]} | wc -w`
        nbf=`expr ${nbf} - 1`

        echo " *** Number of files for Stream ${st} = ${nbf}"


        if [ "${nbf}" != "" ]; then

            jf=0
            while [ ${jf} -lt ${nbf} ]; do
                echo

                jfp1=`expr ${jf} + 1`

                fQ1=${list[${jf}]}
                fQ2=${list[${jfp1}]}

                echo " ==> will use:"; echo "     * ${fQ1}"; echo "     * ${fQ2}"

                ijob=$((ijob+1))

                flog=`basename ${fQ1}`
                flog=`echo ${flog} | sed -e s/".npz"/".out"/g`
                
                CMD="${EXE} ${fQ1} ${fQ2} $((DT_BINS_H*3600/2))"
                echo "  ==> ${CMD}"; echo
                ${CMD} > ${flog} &
                echo; echo; echo
                
                if [ $((ijob%NJPAR)) -eq 0 ]; then
                    echo "Waiting! (ijob = ${ijob})...."
                    wait
                    echo; echo
                fi
                
                jf=`expr ${jf} + 1`

            done

        fi


    done

done
