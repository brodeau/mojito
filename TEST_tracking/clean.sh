#!/bin/bash

. ./conf.bash

if [ "$1" = "1" ]; then
    rm -f figs/tracking/*${NEMO_EXP}*
    rm -f ./nc/*${NEMO_EXP}*_tracking_${YEAR}*.nc npz/Initialized_buoys*.npz

elif [ "$1" = "2" ]; then
    rm -f ./figs/quadgener/*${NEMO_EXP}* npz/[TQ]-mesh_*${NEMO_EXP}*.npz logs/out_S*__*.out

elif [ "$1" = "3" ]; then
    rm -f figs/deformation/*${NEMO_EXP}* npz/DEFORMATIONS_*${NEMO_EXP}* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out
    
elif [ "$1" = "4" ]; then
    rm -f npz/PDF_*${NEMO_EXP}*.npz figs/*PDF*${NEMO_EXP}*.*
    
elif [ "$1" = "all" ]; then
    rm -rf figs npz logs *.png *.nc *.svg *.pdfa nc

else
    echo "Tell me something!"

fi
