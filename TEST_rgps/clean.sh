#!/bin/bash

. ./conf.bash

if [ "$1" = "0" ]; then
    rm -f ./npz/RGPS_batch_selection_dt${DT_BINS_H}h_*${YEAR}*.npz
    
elif [ "$1" = "1" ]; then
    rm -f figs/SELECTION/* nc/SELECTION_*${YEAR}*.nc npz/SELECTION_*${YEAR}*.npz

elif [ "$1" = "c" ]; then
    rm -f figs/coarsify/* nc/SELECTION_*${YEAR}*km.nc

elif [ "$1" = "2" ]; then
    rm -f ./figs/quadgener/* npz/[TQ]-mesh_*.npz logs/*_SELECTION_RGPS_S*.out logs/*_SELECTION_RGPS_S*.err

elif [ "$1" = "3" ]; then
    rm -f figs/deformation/* npz/DEFORMATIONS_* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out
    
elif [ "$1" = "4" ]; then
    rm -f npz/PDF_*.npz figs/*PDF_*.*
    
elif [ "$1" = "all" ]; then
    rm -f *.out *.err *~ \#*
    rm -rf figs logs nc

else
    echo "Tell me something!"

fi
