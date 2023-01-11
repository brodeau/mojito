#!/bin/bash

. ./conf.bash

#EXE="${MOJITO_DIR}/rgps/03_generate_quad_mesh.py"
EXE="${MOJITO_DIR}/rgps/03_generate_quad_mesh_multi.py"

#list=`\ls npz/SELECTION_buoys_RGPS_stream*.npz`
#list=`\ls npz/SELECTION_buoys_RGPS_stream000_*.npz`

pref="npz/SELECTION_buoys_RGPS_stream000"

#for ff in ${list}; do

for dd in ${LIST_RES}; do

    CMD="${EXE} ${pref} ${dd}"
    echo; echo; echo "  ==> ${CMD}"; echo
    ${CMD}
    echo

done

#done
