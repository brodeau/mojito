#!/bin/bash

. ./conf.bash

XTRASFX=""
#XTRASFX="_postQG"


EXE="python3 -u ${SITRCK_DIR}/si3_part_tracker.py"

RESKM="${RESKM}km"

fullStrResKM=""
if [ "${LIST_RD_SS}" = "" ]; then
    fullStrResKM="${RESKM}"
else
    cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
    fullStrResKM="${cr1}-${RESKM}"
fi

usFullRKM="_${fullStrResKM}"

echo
echo " RESKM = ${RESKM}"
echo " fullStrResKM = ${fullStrResKM}"
echo " usFullRKM = ${usFullRKM}"
echo



# 1/ populate the proper NC files to seed from:

if [ "${ISEED_BASE}" = "selection" ]; then
    echo " * Will get RGPS seeding info in: ${DIRIN_PREPARED_RGPS}/nc for RESKM = ${RESKM}"
    list_seed_nc=`\ls ${DIRIN_PREPARED_RGPS}/nc/SELECTION_RGPS_S???_dt${DT_BINS_H}_199?????h??_199?????h??${usFullRKM}${XTRASFX}.nc`

elif [ "${ISEED_BASE}" = "quads" ]; then
    echo " * Will get RGPS seeding info in: ./nc for RESKM = ${RESKM}"
    list_seed_nc=`\ls ./nc/Q-mesh_RGPS_S???_dt${DT_BINS_H}_199?????h??_199?????h??${usFullRKM}${XTRASFX}.nc`

elif [ "${ISEED_BASE}" = "defs" ]; then
    echo " * Will get RGPS seeding info in: ${DIRIN_PREPARED_RGPS}/nc for RESKM = ${RESKM}"
    list_seed_nc=`\ls ${DIRIN_PREPARED_RGPS}/nc/PointsOfQuadsOfDEF_RGPS_S???_dt${DT_BINS_H}_199?????_199?????${usFullRKM}${XTRASFX}.nc`

else
    echo " Unknown value for ISEED_BASE: ${ISEED_BASE}"
    exit
fi



nbf=`echo ${list_seed_nc} | wc -w`

echo " *** We have ${nbf} seeding files !"
echo ${list_seed_nc}
echo

mkdir -p ./logs

ijob=0

for NEMO_EXP in ${LIST_NEMO_EXP}; do

    DIR_FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_EXP}"
    FSI3IN="${DIR_FSI3IN}/${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${SI3DATE1}_${SI3DATE2}_icemod${XTRA_SFX_SI3}.nc4"

    for fnc in ${list_seed_nc}; do

        echo; echo " *** Based on file ${fnc} !"
        SB=`basename ${fnc} | cut -d'_' -f 3`
        if [ "`echo ${SB} | cut -c1-1`" != "S" ]; then
            echo " ERROR: SB does not seem to get the batch string!"; exit
        fi
        lfb=`\ls ${DIRIN_PREPARED_RGPS}/nc/PointsOfQuadsOfDEF_RGPS_${SB}_dt${DT_BINS_H}_199?*_199?*${RESKM}${XTRASFX}.nc`
        nfb=`echo ${lfb} | wc -w`
        echo "      ==> batch is: ${SB}, the ${nfb} files to use:"        
        for fx in ${lfb}; do echo "   -- ${fx}"; done
        echo
        
        if [ "${LIST_RD_SS}" = "" ]; then
            # -i FSI3 -m FMMM -s FSDG [-k KREC] [-e DEND]
            CMD="${EXE} -i ${FSI3IN} -m ${FNMM} -s ${fnc} -N ${NEMO_CONF} -R ${RESKM}" ; # with nc file for init seed...
            echo; echo " *** About to launch:"; echo "     ${CMD}"; echo
            #exit;#lolo
            flog="tracking__${NEMO_EXP}`basename ${fnc} | sed -e s/"SELECTION_RGPS_"/"${NEMO_EXP}_"/g -e s/".nc"/""/g`"
            ${CMD} 1>./logs/${flog}.out 2>./logs/${flog}.err &
            ijob=$((ijob+1))
            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait; echo; echo
            fi
            
        else
            
            for rdss in ${LIST_RD_SS}; do
                if [ "${ISEED_BASE}" = "defs" ]; then
                    fnd=`\ls ${DIRIN_PREPARED_RGPS}/nc/PointsOfQuadsOfDEF_RGPS_${SB}_dt${DT_BINS_H}_199?????_199?????_${rdss}-${RESKM}${XTRASFX}.nc 2>/dev/null`
                    if [ `echo ${fnd} | wc -w` -gt 1 ]; then echo "ERROR: more than 1 candidate for fnd !!!"; echo "  ==> ${fnd}"; exit; fi
                    #echo "LOLO: fnd = ${fnd}"; exit
                else
                    fnd=`echo ${fnc} | sed -e "s|_${cr1}-${RESKM}|_${rdss}-${RESKM}|g"`
                fi
                #
                if [ "${fnd}" != "" ] && [ -f ${fnd} ]; then
                    CMD="${EXE} -i ${FSI3IN} -m ${FNMM} -s ${fnd} -N ${NEMO_CONF} -R ${rdss}-${RESKM}" ; # with nc file for init seed...
                    echo; echo " *** About to launch:"; echo "     ${CMD}"; echo
                    #exit;#lolo
                    flog="tracking__${NEMO_EXP}_`basename ${fnd} | sed -e s/"SELECTION_RGPS_"/"${NEMO_EXP}_"/g -e s/".nc"/""/g`"
                    ${CMD} 1>./logs/${flog}.out 2>./logs/${flog}.err &
                    ijob=$((ijob+1))
                    if [ $((ijob%NJPAR)) -eq 0 ]; then
                        echo "Waiting! (ijob = ${ijob})...."
                        wait; echo; echo
                    fi
                else
                    echo; echo "WARNING: could not find file for res ${rdss}-${RESKM} of batch ${SB} !"; echo
                fi
            done
        fi

    done

done

wait
