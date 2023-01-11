#!/bin/bash

RGPS_DIR="/MEDIA/data/data/RGPS_Kwok_98"

EXE="${HOME}/DEV/mojito/rgps/03_generate_quad_mesh.py"


list=`\ls npz/SELECTION_buoys_RGPS_stream*.npz`

for ff in ${list}; do
    
    for dd in 10 15 20; do
    
        CMD="${EXE} ${ff} ${dd}"
        echo; echo; echo "  ==> ${CMD}"; echo
        ${CMD}
        echo
    done

done

