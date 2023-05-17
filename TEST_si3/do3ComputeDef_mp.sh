#!/bin/bash

#. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

mkdir -p ./logs

ijob=0

for RESKM in ${LCOARSEN[*]}; do

    echo; echo
    echo " *** ${RESKM} km ***"
    echo

    csf="_${RESKM}km"


    for NEMO_EXP in ${LIST_NEMO_EXP}; do
        echo; echo

        # Populating the batches available:
        listQ=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemoTsi3_NoBin_${YEAR}????-??h??t0_${YEAR}????-??h??${csf}.npz`
        #          npz/Q-mesh_NEMO-SI3_NANUK4_BBM00_nemoTsi3_NoBin_19970101t0_19970101_640km.npz
        echo "${listQ}"

        
        list_date_ref=""
        for ff in ${listQ}; do
            date_ref=`echo ${ff} | cut -d_ -f7`
            #echo "LOLO: date_ref = ${date_ref} !"
            list_date_ref+=" ${date_ref}"
        done
        #list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !
        list_date_ref=$(echo "${list_date_ref}" | tr ' ' '\n'| sort -u ) ; # unique only!

        echo; echo " *** List of reference dates:"; echo "${list_date_ref}"; echo
        
        for dr in ${list_date_ref}; do

            echo

            lst=( `\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemoTsi3_NoBin_${dr}_${YEAR}????-??h??${csf}.npz 2>/dev/null` )
            if [ "${lst}" != "" ]; then
                nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "
                #
                if [ ${nf} -eq 2 ]; then
                    fQ1=${lst[0]}
                    fQ2=${lst[1]}
                    echo " ==> will use:"; echo " * ${fQ1}"; echo " * ${fQ2}"
                    flog=`basename ${fQ1}`; flog=`echo ${flog} | sed -e s/".npz"/""/g`; flog="def_${flog}"
                    ijob=$((ijob+1))
                    CMD="${EXE} ${fQ1} ${fQ2} 0 ${MODE}"
                    echo "  ==> ${CMD}"; echo

                    ${CMD} 1>logs/${flog}.out 2>logs/${flog}.err &
                    echo; echo
                    if [ $((ijob%NJPAR)) -eq 0 ]; then
                        echo "Waiting! (ijob = ${ijob})...."
                        wait; echo; echo
                    fi
                fi
            fi


        done

        #done

    done

done

wait
