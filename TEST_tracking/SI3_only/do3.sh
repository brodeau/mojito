#!/bin/bash

# Deformation only
. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

mkdir -p logs


# Populating the batches available:
listQ=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemo_${YEAR}????t0_${YEAR}????_*km.npz`

echo "${listQ}"




echo; echo
echo " *** ${RESKM} km ***"
echo

#  Q-mesh_RGPS_S000_19970104t0_19970104.npz
list=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemo_${YEAR}????t0_${YEAR}????_${RESKM}km.npz`
nbf=`echo ${list} | wc -w`

echo " *** Number of files = ${nbf}"

list_date_ref=""
for ff in ${list}; do
    date_ref=`echo ${ff} | cut -d_ -f6`
    list_date_ref+=" ${date_ref}"
done
list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

echo; echo " *** List of reference dates:"; echo "${list_date_ref}"; echo

for dr in ${list_date_ref}; do
    echo
    lst=( `\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemo_${dr}_${YEAR}????_${RESKM}km.npz` )
    nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "

    if [ ${nf} -eq 2 ]; then

        fQ1=${lst[0]}
        fQ2=${lst[1]}
        echo " ==> will use:"
        echo " * ${fQ1}"
        echo " * ${fQ2}"

        flog=`basename ${fQ1}`
        flog=`echo ${flog} | sed -e s/".npz"/""/g`

        CMD="${EXE} ${fQ1} ${fQ2} 0"
        echo "  ==> ${CMD}"; echo
        ${CMD}
        echo; echo

    fi

done
