#!/bin/bash

. ./conf.bash




for dd in ${LIST_RES}; do

    echo; echo
    echo " *** ${dd} km ***"

    for st in ${LIST_STREAM}; do

        list=( `\ls npz/Q-mesh_stream${st}_*_${dd}km_*.npz` )
        
        echo " ==> will use:"; echo "     * ${list[0]}"; echo "     * ${list[1]}"
        
        ../deformation.py ${list[0]} ${list[1]} 800
        
    done
    
done

