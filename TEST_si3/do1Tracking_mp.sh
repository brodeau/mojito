#!/bin/bash

#. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash


XTRASFX=""
#XTRASFX="_postQG"


EXE="python3 -u ${SITRCK_DIR}/si3_part_tracker.py"


mkdir -p ./logs

ijob=0
    

for RESKM in ${LCOARSEN}; do

    echo; echo; echo ; echo " *** COARSENING: ${RESKM}km !!!"; echo

    
    # 1/ populate the proper NC files to seed from:
    echo " * Will get RGPS seeding info in: ${DIRIN_PREPARED_RGPS} for RESKM = ${RESKM}"
    cxtraRES="_${RESKM}km"

    list_seed_nc=`\ls ./nc/sitrack_seeding_nemoTsi3_${YEAR}*${cxtraRES}.nc`

    nbf=`echo ${list_seed_nc} | wc -w`

    echo " *** We have ${nbf} seeding files !"
    echo ${list_seed_nc}
    echo

    list_date_seed=""
    for ff in ${list_seed_nc}; do
        fb=`basename ${ff}`
        ca=`echo ${fb} | cut -d '.' -f1 | cut -d '_' -f4`
        echo " $fb => $ca"
        list_date_seed+="${ca} "
    done

    echo "List of dates: ${list_date_seed}"

    list_date_seed=( ${list_date_seed} )

    for NEMO_EXP in ${LIST_NEMO_EXP}; do

        DIR_FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_EXP}"
        FSI3IN="${DIR_FSI3IN}/${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${SI3DATE1}_${SI3DATE2}_icemod${XTRA_SFX_SI3}.nc4"

        ii=0
        for cdate in ${list_date_seed[*]}; do

            if [ ${ii} -lt ${nbf} ]; then

                DATE_STOP=${list_date_seed[$((ii+1))]}

                FSEED="./nc/sitrack_seeding_nemoTsi3_${cdate}${cxtraRES}.nc"

                CMD="${EXE} -i ${FSI3IN} -m ${FNMM} -s ${FSEED} -e ${DATE_STOP}"

                echo; echo " *** About to launch:"; echo "     ${CMD}"; echo

                flog="track_${cdate}-${DATE_STOP}_${NEMO_EXP}_${RESKM}"

                ${CMD} 1>./logs/${flog}.out 2>./logs/${flog}.err &
                ijob=$((ijob+1))
                if [ $((ijob%NJPAR)) -eq 0 ]; then
                    echo "Waiting! (ijob = ${ijob})...."
                    wait; echo; echo
                fi

            fi

            ii=$((ii+1))
        done

    done

done
wait