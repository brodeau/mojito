#!/bin/bash

# Deformation only

. ./conf.bash

EXE="${MOJITO_DIR}/deformation.py"

for dd in ${LIST_RES}; do

    echo; echo
    echo " *** ${dd} km ***"
    echo

    for st in ${LIST_STREAM}; do

        list=( `\ls npz/Q-mesh_S${st}_*_${dd}km_*.npz` )
        
        echo " ==> will use:"; echo "     * ${list[0]}"; echo "     * ${list[1]}"
        
        CMD="../deformation.py ${list[0]} ${list[1]} $((DT_BINS_H*3600/2)) ${SIZE_MOSAIC}"
        echo "  ==> ${CMD}"; echo
        ${CMD}
        echo; echo; echo
        
    done
    
done

