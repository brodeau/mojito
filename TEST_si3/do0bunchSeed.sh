#!/bin/bash

. ./conf.bash

EXE="${SITRCK_DIR}/tools/generate_idealized_seeding.py"

EXE2="${SITRCK_DIR}/tools/listofdates.py"

LIST_DATES=`${EXE2} ${DATE1} ${DATE2} ${DT_INC_DAYS}`

echo
echo "LIST_DATES = ${LIST_DATES}"
echo

NEMO_EXP=`echo ${LIST_NEMO_EXP} | cut -d' ' -f1`

DIR_FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_EXP}"
FSI3IN="${DIR_FSI3IN}/${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${SI3DATE1}_${SI3DATE2}_icemod${XTRA_SFX_SI3}.nc4"

mkdir -p ./logs

ijob=0


for date in ${LIST_DATES}; do

    YYYY=`echo ${date} | cut -c1-4`
    MM=`echo ${date} | cut -c5-6`
    DD=`echo ${date} | cut -c7-8`
    hh=`echo ${date} | cut -c10-11`


    #echo " DD = ${DD} / hh = ${hh}"; exit
    
    export NDATE0="${YYYY}-${MM}-${DD}_${hh}"
    export LDATE0="${NDATE0}:00:00"


    echo; echo ; echo ; echo " #### DATE: $NDATE0 ####"; echo
    
    for RESKM in ${LCOARSEN}; do

        echo; echo; echo " *** RESKM = ${RESKM} ***"        

        str="${NDATE0}_${RESKM}km"

        flog="seeding_${str}"

        fout="./nc/sitrack_seeding_nemoTsi3_${NDATE0}_${RESKM}km.nc"

        if [ ! -f ${fout} ]; then

            CMD="${EXE} -d ${LDATE0} -m ${FNMM} -i ${FSI3IN} -k 0 -f ${FFSM} -C ${RESKM}"
            echo
            echo " *** About to launch:"; echo "     ${CMD}"; echo
            ijob=$((ijob+1))
            ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
            echo
            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait; echo; echo
            fi

        fi
    done
done
