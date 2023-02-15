#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/generate_quad_mesh.py"

if [ "$1" = "" ]; then
    echo "USAGE: $0 <file_pos_mojito.nc>"
    exit
fi

FIN="$1"

# Number of records inside netCDF file:
Nr0=`ncdump -h ${FIN} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
echo
echo " * ${Nr0} records in ${FIN}!"
Nr=$(((Nr0+1)/DT_BINS_H))
Nr=$((Nr+1))
echo "   => subsampling at ${DT_BINS_H} hours give ${Nr} useful records!"
echo $((Nr*DT_BINS_H))
#exit


ijob=0

#CMD="${EXE} ${FIN} 0,1 10" ; # RGPS seed

RESKM=10

mkdir -p logs

for ii in $(seq 0 $((Nr-2))); do
    echo $ii
    #
    rec1=$((ii*DT_BINS_H))
    rec2=$(((ii+1)*DT_BINS_H))

    echo " * rec1, rec2 = ${rec1} ${rec2}"

    cflog="logs/out_${ii}_$((ii+1)).out"

    #if [ ! -f ${cfQ1} ] || [ ! -f ${cfQ2} ]; then
        ijob=$((ijob+1))
        echo " *** Construction of Quadrangles"
        CMD="${EXE} ${FIN} ${rec1},${rec2} ${RESKM}" ; # RGPS seed
        echo "  ==> ${CMD}"; echo
        ${CMD} > ${cflog} &
    echo; sleep 1 ; echo
    #else
    #    echo; echo
    #    echo " Skipping generation of Quads as ${cfQ1} & ${cfQ2} already there!"
    #    echo; echo
    #fi

    if [ $((ijob%NJPAR)) -eq 0 ]; then
        echo "Waiting! (ijob = ${ijob})...."
        wait
        echo; echo
    fi

done

wait

echo
echo " *** `date` ALL done!"
echo
