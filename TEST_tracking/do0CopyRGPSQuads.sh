#!/bin/bash

. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
. ./conf.bash

EXE="python3 -u /home/laurent/DEV/mojito/util/convert_Qnpz2SeedNC.py"

ijob=0
mkdir -p logs

if [ "${ISEED_BASE}" != "quads" ]; then
    echo "There is no need to call me since ISEED_BASE is not 'quads' !!!"
fi



# 1/ populate the proper npz Quad files to seed from:
echo " * Will get RGPS seeding info in: ${DIRIN_PREPARED_RGPS}/npz for RESKM = ${RESKM}"
cxtraRES=""
if [ "${LIST_RD_SS}" = "" ]; then
    cxtraRES="_${RESKM}km"
else
    cr1=`echo ${LIST_RD_SS} | cut -d' ' -f1` ; # premiere resolution `rd_ss` !!!
    cxtraRES="_${cr1}-${RESKM}km"
fi


echo " RESKM = ${RESKM}"
echo " cxtraRES = ${cxtraRES}"

#list_seed_qnpz=`\ls ${DIRIN_PREPARED_RGPS}/SELECTION_RGPS_S???_dt${DT_BINS_H}_${YEAR}????h??_${YEAR}????h??${cxtraRES}${XTRASFX}.nc`
list_seed_qnpz=`\ls ${DIRIN_PREPARED_RGPS}/npz/Q-mesh_RGPS_S???_dt${DT_BINS_H}_${YEAR}????-??h??t0_${YEAR}????-??h??${cxtraRES}${XTRASFX}.npz`

# need to keep only the first occurence of each:
lnew=""
icpt=0
for ff in ${list_seed_qnpz}; do
    if [ $((icpt%2)) -eq 0 ]; then
        lnew+="${ff} "
    fi
    icpt=$((icpt+1))
done
list_seed_qnpz=${lnew}
nbf=`echo ${list_seed_qnpz} | wc -w`

echo " *** We have ${nbf} seeding files !"
echo ${list_seed_qnpz}
echo


for ff in ${list_seed_qnpz}; do

    echo " *** Will use ${ff} !"
    fd=`dirname ${ff}`
    fb=`basename ${ff}`
    croot=`echo ${fb} | cut -d'-' -f1-2`
    echo ${croot}

    flog="0getfromQuads_${croot}${cxtraRES}${XTRASFX}"

    ijob=$((ijob+1))
    CMD="${EXE} ${fd}/${croot}"
    echo; echo " *** About to launch:"; echo "     ${CMD}"; echo

    ${CMD} 1>"./logs/${flog}.out" 2>"./logs/${flog}.err" &

    if [ $((ijob%NJPAR)) -eq 0 ]; then
        echo "Waiting! (ijob = ${ijob})...."
        wait
        echo; echo
    fi

    icpt=$((icpt+1))

done

exit
rm -f *.png
