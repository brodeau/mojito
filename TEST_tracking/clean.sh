#!/bin/bash

. ./conf.bash



st1="_${RESKM}km"
st2="-${RESKM}km"

if [ "$1" = "1" ]; then
    rm -rf figs/tracking/${RESKM}km
    rm -f ./nc/NEMO-SI3_*_tracking*${st1}*.nc
    rm -f npz/Initialized_buoys_PointsOfQuadsOfDEF_*${st1}_*.npz npz/Initialized_buoys_PointsOfQuadsOfDEF_*${st2}_*.npz
    rm -f ./logs/tracking__*${st1}*
    rm -f ./nc/NEMO-SI3_*_tracking*${st2}*.nc
    rm -f ./logs/tracking__*${st2}*

elif [ "$1" = "2" ]; then
    rm -rf ./figs/quadgener/${RESKM}km
    rm -f npz/[TQ]-mesh_*${st1}*.npz logs/quadgener_*${st1}*.*
    rm -f npz/[TQ]-mesh_*${st2}*.npz logs/quadgener_*${st2}*.*

elif [ "$1" = "3" ]; then
    rm -rf figs/deformation/${RESKM}km
    rm -f npz/DEFORMATIONS_NEMO-SI3*${st1}* npz/DEFORMATIONS_NEMO-SI3*${st2}*
    rm -f npz/DEFORMATIONS_NEMO-SI3*_SLCT${RESKM}km.npz
    rm -f logs/err_Q-mesh_*${st1}*.err logs/out_Q-mesh_*${st1}*.out logs/def__*${st1}*
    rm -f logs/err_Q-mesh_*${st2}*.err logs/out_Q-mesh_*${st2}*.out logs/def__*${st2}*
    
elif [ "$1" = "4" ]; then
    rm -f npz/def_*${st1}*  logs/def__*${st1}*
    rm -f npz/def_*${st2}*  logs/def__*${st2}*
    
elif [ "$1" = "5" ]; then
    rm -f npz/PDF_*${st1}*.npz figs/*PDF_*${st1}*.*
    rm -f npz/PDF_*${st2}*.npz figs/*PDF_*${st2}*.*
    
elif [ "$1" = "all" ]; then
    rm -f *.out *.err *~ \#*
    rm -rf figs logs nc
    rm -f npz/PDF_* npz/DEFORMATIONS_* npz/[QT]-mesh* npz/def_*

else
    echo "Tell me something!"

fi
