#!/bin/bash

. ./conf.bash


if [ "$1" = "" ]; then
    echo "USAGE: $0 <Nb. batches>"
    exit
fi

Nb=${1}

#dir1=`pwd`
dir1="." ; cnm1="NEMO-SI3_NANUK4_BBM2302"
dir2="../TEST_rgps" ; cnm2="RGPS"



jb=0

while [ ${jb} -le ${Nb} ]; do

    sb="S`printf "%03d" ${jb}`"

    echo
    echo " *** Looking at batch # ${sb}"


    for rskm in ${LIST_RD_SS}; do

        echo "  ==> ${rskm} km"

        fdef1=`\ls ${dir1}/npz/DEFORMATIONS_${cnm1}_${sb}_dt72_1997*_${rskm}-${RESKM}km.npz 2>/dev/null`
        nbf1=`echo ${fdef1} | wc -w`
        echo "    + File 1: $fdef1 ; nl=  $nbf1"
        
        fdef2=`\ls ${dir2}/npz/DEFORMATIONS_${cnm2}_${sb}_dt72_1997*_${rskm}-${RESKM}km.npz 2>/dev/null`
        nbf2=`echo ${fdef2} | wc -w`
        echo "    + File 2: $fdef2 ; nl=  $nbf2"

        if [ "${fdef1}" = "" ] && [ ${nbf2} -eq 1 ] && [ -f ${fdef2} ]; then
            echo "PROBLEM: file for ${cnm1} is missing while that of ${cnm2} is here!!!"; exit
        fi


        
        if [ ${nbf1} -ne ${nbf2} ]; then
            echo "PROBLEM: not same files!!!"; exit
        fi



        

        
        echo; echo
    done
    


    jb=`expr ${jb} + 1`
done
