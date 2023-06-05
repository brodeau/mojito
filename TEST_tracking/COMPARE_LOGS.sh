#!/bin/bash

. ./conf.bash


if [ "$2" = "" ]; then
    echo "USAGE: $0 <first> <Nb. batches>"
    exit
fi

istart=${1}
Nb=${2}

#dir1=`pwd`
dir1="." ; cnm1="NEMO-SI3_NANUK4_BBM2302"
dir2="../TEST_rgps" ; cnm2="RGPS"



jb=${istart}

while [ ${jb} -le ${Nb} ]; do

    sb="S`printf "%03d" ${jb}`"

    echo
    echo " *** Looking at batch # ${sb}"


    for rskm in ${LIST_RD_SS}; do

        echo "  ==> BATCH ${sb} / ${rskm}-${RESKM}km"

        fdef1=`\ls ${dir1}/npz/DEFORMATIONS_${cnm1}_${sb}_dt72_1997*_${rskm}-${RESKM}km.npz 2>/dev/null`
        nbf1=`echo ${fdef1} | wc -w`
        echo "    + File 1: $fdef1 ; nl=  $nbf1"
        
        fdef2=`\ls ${dir2}/npz/DEFORMATIONS_${cnm2}_${sb}_dt72_1997*_${rskm}-${RESKM}km.npz 2>/dev/null`
        nbf2=`echo ${fdef2} | wc -w`
        echo "    + File 2: $fdef2 ; nl=  $nbf2"

        if [ "${fdef1}" = "" ] && [ ${nbf2} -eq 1 ] && [ -f ${fdef2} ]; then
            echo "PROBLEM: file for ${cnm1} is missing while that of ${cnm2} is here!!!"; exit
        fi

        if [ "${fdef2}" = "" ] && [ ${nbf1} -eq 1 ] && [ -f ${fdef1} ]; then
            echo "PROBLEM: file for ${cnm2} is missing while that of ${cnm1} is here!!!"; exit
        fi
        
        if [ ${nbf1} -ne ${nbf2} ]; then
            echo "PROBLEM: not same files!!!"; exit
        fi



        if [ ${nbf1} -gt 1 ] || [ ${nbf2} -gt 1 ]; then
            echo "ERROR: too many files: ${nbf1} , ${nbf2}"; exit
        fi


        if [ ${nbf1} -eq 1 ]; then
            nP1=`${MOJITO_DIR}/util/get_number_of_defs.py ${fdef1}`
            nP2=`${MOJITO_DIR}/util/get_number_of_defs.py ${fdef2}`


            echo " Nb. of deformations in each file: ${nP1}, ${nP2}"

            if [ ${nP1} -ne ${nP2} ]; then
                echo "PROBLEM: not same number of deformtion points!!!"; exit
            fi

        fi
        

        
        echo; echo
    done
    


    jb=`expr ${jb} + 1`
done
