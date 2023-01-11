#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/rgps/03_generate_quad_mesh.py"

list=`\ls npz/SELECTION_buoys_RGPS_stream*.npz`

for ff in ${list}; do
    
    for dd in ${LIST_RES}; do
    
        CMD="${EXE} ${ff} ${dd}"
        echo; echo; echo "  ==> ${CMD}"; echo
        ${CMD}
        echo
        
    done

done

