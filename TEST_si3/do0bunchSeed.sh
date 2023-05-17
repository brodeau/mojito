#!/bin/bash

. ./conf.bash

EXE="${SITRCK_DIR}/tools/generate_idealized_seeding.py"

EXE2="${SITRCK_DIR}/tools/listofdates.py"

NEMO_EXP=`echo ${LIST_NEMO_EXP} | cut -d' ' -f1`

DIR_FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_EXP}"
FSI3IN="${DIR_FSI3IN}/${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${SI3DATE1}_${SI3DATE2}_icemod${XTRA_SFX_SI3}.nc4"

mkdir -p ./logs

ijob=0

ires=-1

for RESKM in ${LCOARSEN[*]}; do
    ires=$((ires+1))
    DT_INC_DAYS=${LDTINCRM[${ires}]}

    echo; echo; echo " *** RESKM = ${RESKM} => DT = ${DT_INC_DAYS} ***"

    LIST_DATES=`${EXE2} ${DATE1} ${DATE2} ${DT_INC_DAYS}`

    echo "    ===> LIST_DATES = ${LIST_DATES}"; echo

    idt=0
    for date in ${LIST_DATES}; do

        idt=`expr ${idt} + 1`

        YYYY=`echo ${date} | cut -c1-4`
        MM=`echo ${date} | cut -c5-6`
        DD=`echo ${date} | cut -c7-8`
        hh=`echo ${date} | cut -c10-11`

        export NDATE0="${YYYY}-${MM}-${DD}_${hh}"
        export LDATE0="${NDATE0}:00:00"

        echo; echo ; echo ; echo " #### DATE: $NDATE0 ####"; echo

        str="${NDATE0}_${RESKM}km"

        flog="seeding_${str}"

        fmaskRGPS=${FFSM}
        if [ "${DT_INC_DAYS}" = "0.25" ]; then
            case ${idt} in
                2)
                    fmaskRGPS=`echo ${FFSM} | sed -e "s|_all.nc|_AlskGrnlnd.nc|g"` ;;
                3)
                    fmaskRGPS=`echo ${FFSM} | sed -e "s|_all.nc|_GrnlndNwSbr.nc|g"` ;;
                4)
                    fmaskRGPS=`echo ${FFSM} | sed -e "s|_all.nc|_CndNns.nc|g"`
                    idt=0 ;;
            esac
        fi

        sdate0=`echo ${NDATE0} | sed -e s/'-'/''/g`
        fout="./nc/sitrack_seeding_nemoTsi3_${sdate0}_${RESKM}km.nc"

        if [ ! -f ${fout} ]; then

            CMD="${EXE} -d ${LDATE0} -m ${FNMM} -i ${FSI3IN} -k 0 -f ${fmaskRGPS} -C ${RESKM}"
            echo
            echo " *** About to launch:"; echo "     ${CMD}"; echo
            ijob=$((ijob+1))
            ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &
            echo
            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait; echo; echo
            fi

        else
            echo " File ${fout} is already there! Nothing to do!"
        fi

    done

done
