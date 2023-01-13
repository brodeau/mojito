#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/rgps/03_generate_quad_mesh_multi.py"


for st in ${LIST_STREAM}; do


    vlist=( `\ls npz/SELECTION_buoys_RGPS_stream${st}*.npz` )

    #echo "${vlist[*]}"; echo
    nf=`echo "${vlist[*]}" | wc -w`
    echo; echo "*** ${nf} files!"; echo


    for ii in $(seq 0 $((nf-2))); do
        echo $ii
        fref=${vlist[${ii}]}
        ftst=${vlist[$((ii+1))]}

        echo " => work with: ${fref} & ${ftst}"

        CMD="${EXE} ${fref},${ftst} ${dd}"
        echo "  ==> ${CMD}"; echo
        ${CMD}
        echo; echo
    done


done
exit





for st in ${LIST_STREAM}; do

    pref="npz/SELECTION_buoys_RGPS_stream${st}"


    for dd in ${LIST_RES}; do


    done

done
