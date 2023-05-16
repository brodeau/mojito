#!/bin/bash

. ./conf.bash

if [ "$1" = "1" ]; then
    rm -rf figs/tracking
    for NEMO_EXP in ${LIST_NEMO_EXP}; do
        rm -f ./nc/NEMO-SI3_NANUK4_${NEMO_EXP}_tracking_nemoTsi3_*.nc
    done

elif [ "$1" = "2" ]; then
    for NEMO_EXP in ${LIST_NEMO_EXP}; do
        rm -f ./nc/NEMO-SI3_NANUK4_${NEMO_EXP}_tracking_nemoTsi3_idlSeed_*km.nc
    done
    rm -f ./figs/quadgener/* npz/[TQ]-mesh_*.npz logs/quadgener_*

elif [ "$1" = "3" ]; then
    rm -rf figs/deformation/* npz/DEFORMATIONS_* logs/def_*
    
elif [ "$1" = "4" ]; then
    rm -f npz/def_*
    
elif [ "$1" = "5" ]; then
    rm -f npz/PDF_*.npz figs/*PDF_*.*
    
elif [ "$1" = "all" ]; then
    rm -f *.out *.err *~ \#*
    rm -rf figs logs nc npz
    
else
    echo "Tell me something!"

fi
