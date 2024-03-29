#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/generate_quad_mesh.py"

if [ "$1" = "" ]; then
    echo "USAGE: $0 <file_pos_mojito.nc>"
    exit
fi

FIN="$1"

# Number of records inside netCDF file:
Nr=`ncdump -h ${FIN} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
echo
echo " * ${Nr} records in ${FIN}!"

ijob=0

RESKM=10

mkdir -p logs

for ii in $(seq 0 $((Nr-2))); do
    echo $ii
    
    cflog="logs/out_${ii}_$((ii+1)).out"
    
    ijob=$((ijob+1))
    echo " *** Construction of Quadrangles"
    CMD="${EXE} ${FIN} ${ii},$((ii+1)) ${RESKM}" ; # RGPS seed
    echo "  ==> ${CMD}"; echo
    ${CMD} > ${cflog} &
    echo; sleep 1 ; echo
    
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
