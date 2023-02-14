#!/bin/bash

. ./conf.bash


EXE="${MOJITO_DIR}/rgps/03_generate_quad_mesh_multi.py"

mkdir -p logs

ijob=0


for dd in ${LIST_RES}; do

    csuff="Sampled"
    #if [ ${dd} -eq 10 ]; then csuff="NoSample"; fi

    istr=${NSTREAMA}
    while [ ${istr} -le ${NSTREAMB} ]; do
              
        echo; echo; echo

        cstr=`printf "%03d" ${istr}`

        vlist=( `\ls npz/SELECTION_buoys_RGPS_S${cstr}*.npz` )

        #echo "${vlist[*]}"; echo
        nf=`echo "${vlist[*]}" | wc -w`
        echo; echo "*** ${nf} files!"; echo


        for ii in $(seq 0 $((nf-2))); do
            echo $ii
            fref=${vlist[${ii}]}      ; cf1=`basename ${fref} | cut -d'.' -f1`
            ftst=${vlist[$((ii+1))]}  ; cf2=`basename ${ftst} | cut -d'.' -f1`

            echo; echo " => work with: ${fref} & ${ftst}"; echo

            cdt1=`echo ${cf1} | cut -d_ -f5`
            cdt2=`echo ${cf2} | cut -d_ -f5`
            chr1=`echo ${cf1} | cut -d_ -f6`
            chr2=`echo ${cf2} | cut -d_ -f6`

            # The two quadrangle files to be generated:
            cfQ1="npz/Q-mesh_S${cstr}_${cdt1}t0_${cdt1}_${chr1}_${dd}km_${csuff}.npz"
            cfQ2="npz/Q-mesh_S${cstr}_${cdt1}t0_${cdt2}_${chr2}_${dd}km_${csuff}.npz"
            #
            cflog="logs/out_S${cstr}_${cdt1}_${chr1}__${cdt2}_${chr2}_${dd}km_${csuff}.out"
            
            if [ ! -f ${cfQ1} ] || [ ! -f ${cfQ2} ]; then
                ijob=$((ijob+1))
                echo " *** Construction of Quadrangles"
                CMD="${EXE} ${fref},${ftst} ${dd}"
                echo "  ==> ${CMD}"; echo
                ${CMD} > ${cflog} &
                echo; echo
            else
                echo; echo
                echo " Skipping generation of Quads as ${cfQ1} & ${cfQ2} already there!"
                echo; echo
            fi

            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait
                echo; echo
            fi
            
        done

        istr=`expr ${istr} + 1`
    done
done

wait

echo
echo " *** `date` ALL done!"
echo