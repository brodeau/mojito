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

#CMD="${EXE} ${FIN} 0,1 10" ; # RGPS seed

RESKM=10

mkdir -p logs

for ii in $(seq 0 $((Nr-2))); do
    echo $ii

    #fref=${vlist[${ii}]}      ; cf1=`basename ${fref} | cut -d'.' -f1`
    #ftst=${vlist[$((ii+1))]}  ; cf2=`basename ${ftst} | cut -d'.' -f1`

    #echo; echo " => work with: ${fref} & ${ftst}"; echo

    #cdt1=`echo ${cf1} | cut -d_ -f5`
    #cdt2=`echo ${cf2} | cut -d_ -f5`
    #chr1=`echo ${cf1} | cut -d_ -f6`
    #chr2=`echo ${cf2} | cut -d_ -f6`

    # The two quadrangle files to be generated:
    #cfQ1="npz/Q-mesh_S${cstr}_${cdt1}t0_${cdt1}_${chr1}_${dd}km_${csuff}.npz"
    #cfQ2="npz/Q-mesh_S${cstr}_${cdt1}t0_${cdt2}_${chr2}_${dd}km_${csuff}.npz"
    #
    cflog="logs/out_${ii}_$((ii+1)).out"

    #if [ ! -f ${cfQ1} ] || [ ! -f ${cfQ2} ]; then
        ijob=$((ijob+1))
        echo " *** Construction of Quadrangles"
        CMD="${EXE} ${FIN} ${ii},$((ii+1)) ${RESKM}" ; # RGPS seed
        echo "  ==> ${CMD}"; echo
        ${CMD} > ${cflog} &
        echo; echo
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


#wait

echo
echo " *** `date` ALL done!"
echo
