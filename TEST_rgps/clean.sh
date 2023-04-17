#!/bin/bash

. ./conf.bash

if [ "$1" = "0" ]; then
    rm -f ./npz/RGPS_batch_selection_dt${DT_BINS_H}h_*${YEAR}*.npz
    
elif [ "$1" = "1" ]; then
    rm -f figs/SELECTION/* nc/SELECTION_*${YEAR}*.nc npz/SELECTION_*${YEAR}*.npz

elif [ "$1" = "c" ]; then
    rm -f figs/coarsify/* nc/SELECTION_*${YEAR}*km.nc

elif [ "$1" = "2" ]; then
    rm -f ./figs/quadgener/* npz/[TQ]-mesh_*.npz logs/*_SELECTION_RGPS_S*.out logs/genquads_*
    rm -f ./nc/*_postQG.nc

elif [ "$1" = "3" ]; then
    rm -f figs/deformation/* npz/DEFORMATIONS_* logs/def_*
    
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
