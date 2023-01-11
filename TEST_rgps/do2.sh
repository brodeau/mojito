#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/rgps/03_generate_quad_mesh_multi.py"

for st in ${LIST_STREAM}; do

    pref="npz/SELECTION_buoys_RGPS_stream${st}"


    for dd in ${LIST_RES}; do
        
        CMD="${EXE} ${pref} ${dd}"
        echo; echo; echo "  ==> ${CMD}"; echo
        ${CMD}
        echo

    done

done
