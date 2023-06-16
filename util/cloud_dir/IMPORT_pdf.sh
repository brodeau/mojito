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



# /home/laurent/tmp/MOJITO/TEST_tracking

HERE=`pwd`

#for cw in "bbm" "evp"; do
for EXP in ${lst_exps}; do



    cw=${EXP,,}


    echo " ${cw}, ${EXP}"

    
    mkdir -p ${cw}
    
    cd ${cw}/

    echo; echo " * Inside `pwd` !"

    csfx=${cw}

    #echo "LOLO: cw = ${cw}"
    
    if [ "${cw}" != "rgps" ]; then
        #csfx="si3"
        csfx="tracking"
    fi
        
    dir=${DIR_TRCK_REMOTE}

    #echo " ${dir}/TEST_${csfx}/${DIR_NPZ}"
    #exit

    echo " * Will import from: ${dir}/TEST_${csfx}/${DIR_NPZ} !"
    
    import_frazilo ${dir}/TEST_${csfx}/${DIR_NPZ}/def_\*${EXP}\*.npz .
    import_frazilo ${dir}/TEST_${csfx}/${DIR_NPZ}/PDF_\*${EXP}\*.npz .
    
    cd ../

done
