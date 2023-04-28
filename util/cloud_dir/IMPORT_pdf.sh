#!/bin/bash

. ./conf.bash

lst_exps=`echo ${EXPS} | sed -e s/','/' '/g -e s/'BBM00'/''/g`
echo " => ${lst_exps} "


import_frazilo()
{
    rsync -avP -e "ssh -p ${PORT_FRAZILO} -c aes128-ctr" laurent@localhost:${1} ${2}
}


DIR_NPZ="npz"
#DIR_NPZ="BAK_DEF"

#if [ "${1}" = "c" ]; then
#    for cw in "rgps" "tracking"; do rm -rf TEST_${cw}/npz/PDF*.npz; done
#    exit
#fi


#DIR_RGPS_REMOTE="/home/laurent/DEV/mojito"
#DIR_TRCK_REMOTE="/home/laurent/tmp/MOJITO"
DIR_TRCK_REMOTE="/home/laurent/DEV/mojito"


# /home/laurent/tmp/MOJITO/TEST_tracking

HERE=`pwd`

#for cw in "bbm" "evp"; do
for EXP in ${lst_exps}; do



    cw=${EXP,,}


    echo " ${cw}, ${EXP}"

    
    mkdir -p ${cw}
    
    cd ${cw}/

    echo; echo " * Inside `pwd` !"

    #dir=${DIR_RGPS_REMOTE}
    #csfx=${cw}
    #cxtr="RGPS"
    #if [ "${cw}" != "rgps" ]; then
    dir=${DIR_TRCK_REMOTE}
    csfx="si3"
    #cxtr="${cw^^}0"
    #cxtr="${cw}"
    #fi
    
    import_frazilo ${dir}/TEST_${csfx}/${DIR_NPZ}/def_\*${EXP}\*.npz .
    import_frazilo ${dir}/TEST_${csfx}/${DIR_NPZ}/PDF_\*${EXP}\*.npz .
    
    cd ../

done
