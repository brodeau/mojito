#!/bin/bash

. ./conf.bash

if [ "$1" = "0" ]; then
    rm -f ./npz/RGPS_batch_selection_dt${DT_BINS_H}h_*${YEAR}*.npz
    rm -rf figs/SELECTION
    rm -f nc/SELECTION_*.nc

elif [ "$1" = "1" ]; then
    rm -rf figs/coarsify
    rm -f nc/SELECTION_*km.nc logs/coarsify_*

elif [ "$1" = "2" ]; then
    rm -rf ./figs/quadgener
    rm -f npz/[TQ]-mesh_*.npz logs/genquads_*
    rm -f ./nc/*_postQG.nc

elif [ "$1" = "3" ]; then
    rm -rf figs/deformation
    rm -f npz/DEFORMATIONS_* logs/def_*
    rm -f nc/PointsOfQuadsOfDEF_*.nc npz/QUADSofDEF_*
    rm -f npz/QUADSofDEF_*
    
elif [ "$1" = "4" ]; then
    rm -f npz/def_*
    
elif [ "$1" = "5" ]; then
    rm -f npz/PDF_*.npz figs/*PDF_*.*
    
elif [ "$1" = "all" ]; then
    rm -f *.out *.err *~ \#*
    rm -rf figs logs nc
    rm -f npz/PDF_* npz/DEFORMATIONS_* npz/[QT]-mesh* npz/def_* npz/QUADSofDEF_*

elif [ "$1" = "less" ]; then
    rm -f *.out *.err *~ \#*
    rm -rf figs logs nc/SELECTION_*km.nc
    rm -f npz/PDF_* npz/DEFORMATIONS_* npz/[QT]-mesh* npz/def_*

    
else
    echo "Tell me something!"

fi
