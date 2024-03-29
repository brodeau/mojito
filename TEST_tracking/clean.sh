#!/bin/bash

. ./conf.bash

if [ "$1" = "1" ]; then
    rm -f figs/tracking/*${NEMO_EXP}*
    rm -f ./nc/NEMO-SI3_NANUK4_${NEMO_EXP}_tracking_S???_dt*_*km.nc

elif [ "$1" = "2" ]; then
    rm -f ./figs/quadgener/* msh/[TQ]-mesh_*.npz logs/*_SELECTION_RGPS_S*.out logs/*_SELECTION_RGPS_S*.err

elif [ "$1" = "3" ]; then
    rm -f figs/deformation/* npz/DEFORMATIONS_* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out
    
elif [ "$1" = "4" ]; then
    rm -f npz/def_*
    
elif [ "$1" = "5" ]; then
    rm -f npz/PDF_*.npz figs/*PDF_*.*
    
elif [ "$1" = "all" ]; then
    rm -f *.out *.err *~ \#* out
    rm -rf figs logs nc npz msh

else
    echo "Tell me something!"

fi
