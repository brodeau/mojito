#!/bin/bash

. ./conf.bash

i_do_def=0 ; # computes deformations or not...

EXE1="${MOJITO_DIR}/rgps/03_generate_quad_mesh_multi.py"
EXE2="${MOJITO_DIR}/deformation.py"


for dd in ${LIST_RES}; do

    csuff="Sampled"
    #if [ ${dd} -eq 10 ]; then csuff="NoSample"; fi

    for st in ${LIST_STREAM}; do


        vlist=( `\ls npz/SELECTION_buoys_RGPS_stream${st}*.npz` )

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
            cfQ1="npz/Q-mesh_S${st}_${cdt1}t0_${cdt1}_${chr1}_${dd}km_${csuff}.npz"
            cfQ2="npz/Q-mesh_S${st}_${cdt1}t0_${cdt2}_${chr2}_${dd}km_${csuff}.npz"
            #
            
            if [ ! -f ${cfQ1} ] || [ ! -f ${cfQ2} ]; then
                echo " *** Construction of Quadrangles"
                CMD="${EXE1} ${fref},${ftst} ${dd}"
                echo "  ==> ${CMD}"; echo
                ${CMD}
                echo; echo
            else
                echo; echo
                echo " Skipping generation of Quads as ${cfQ1} & ${cfQ2} already there!"
                echo; echo
            fi

            if [ ${i_do_def} -eq 1 ]; then

                echo " *** Computing deformations"
                if [ -f ./${cfQ1} ] && [ -f ./${cfQ2} ]; then
                    echo " *** Computation of deformartions"
                    CMD="${EXE2} ./${cfQ1} ./${cfQ2} ${SIZE_MOSAIC}"
                    echo "  ==> ${CMD}"; echo
                    ${CMD}
                    echo; echo
                else
                    echo "PROBLEM: missing Quad files:"
                    echo " ./${cfQ1}  or  ./${cfQ2} "
                    echo
                    exit
                fi

            fi

        done


    done


done
