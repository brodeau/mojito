#!/bin/bash

. ./conf.bash

if [ "$1" = "1" ]; then
    rm -rf figs/tracking/*
    rm -f ./nc/NEMO-SI3_NANUK4_${NEMO_EXP}_tracking_S???_dt*_*km.nc

elif [ "$1" = "2" ]; then
    rm -f ./figs/quadgener/* npz/[TQ]-mesh_*.npz logs/*_SELECTION_RGPS_S*.out logs/*_SELECTION_RGPS_S*.err

elif [ "$1" = "3" ]; then
    rm -rf figs/deformation/*
    rm -f npz/DEFORMATIONS_* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out
    
elif [ "$1" = "4" ]; then
    rm -f npz/def_*
    
elif [ "$1" = "5" ]; then
    rm -f npz/PDF_*.npz figs/*PDF_*.*
    
elif [ "$1" = "all" ]; then
    rm -f *.out *.err *~ \#*
    rm -rf figs logs nc
    rm -f npz/PDF_* npz/DEFORMATIONS_* npz/[QT]-mesh* npz/def_*

else
    echo "Tell me something!"

fi
