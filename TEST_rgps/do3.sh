#!/bin/bash



# npz/Q-mesh_stream001_19970106_18h00_19970106_18h00_10km_Sampled.npz
# npz/Q-mesh_stream001_19970106_18h00_19970109_18h00_10km_Sampled.npz



CSTREAMS="stream000_19970104_18h00 stream001_19970106_18h00_"

CPREF="npz/Q-mesh_stream000_19970104_18h00 npz/Q-mesh_stream001_19970106_18h00_"



for dd in ${LIST_RES}; do

    echo; echo
    echo " *** ${dd} km ***"
    
    for cstr in ${CSTREAMS}; do
        echo
        echo " * Stream ${cstr}"
        
        CPREF="Q-mesh_${cstr}"

        list=( `\ls npz/${CPREF}*_${dd}km_*.npz` )

        
        echo " ==> will use:"; echo "     * ${list[0]}"; echo "     * ${list[1]}"
        
        ../deformation.py ${list[0]} ${list[1]} 120

    done
    
done

