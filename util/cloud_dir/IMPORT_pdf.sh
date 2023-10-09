#!/bin/bash

. ./conf.bash

lst_exps=`echo ${EXPS} | sed -e s/','/' '/g -e s/'BBM00'/''/g`
echo " => ${lst_exps} "


import_frazilo()
{
    rsync -avP -e "ssh -p ${PORT_FRAZILO} -c aes128-ctr" laurent@localhost:${1} ${2}
}


DIR_NPZ="npz"

HERE=`pwd`


lscales=`echo ${SCALES} | sed -e s/','/' '/g`

#for cw in "bbm" "evp"; do
for EXP in ${lst_exps}; do


    for SCL in ${lscales}; do


        cw=${EXP,,}


        echo; echo " ${cw}, ${EXP}, ${SCL}"

        
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


        echo " * Will import from: ${dir}/TEST_${csfx}/${DIR_NPZ} !"

        import_frazilo ${dir}/TEST_${csfx}/${DIR_NPZ}/def_\*${EXP}\*_${SCL}km.npz .
        #import_frazilo ${dir}/TEST_${csfx}/${DIR_NPZ}/PDF_\*${EXP}\*.npz .

        cd ../

    done
done
