#!/bin/bash

. ./conf.bash

if [ "$1" = "0" ]; then
    rm -f ./npz/RGPS_*${YEAR}*.npz
    
elif [ "$1" = "1" ]; then
    rm -f figs/SELECTION/*

elif [ "$1" = "2" ]; then
    rm -f ./figs/quadgener/* npz/[TQ]-mesh_*.npz logs/out_S*__*.out

elif [ "$1" = "3" ]; then
    rm -f figs/deformation/* npz/DEFORMATIONS_* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out
    
elif [ "$1" = "4" ]; then
    rm -f npz/PDF_*.npz *PDF*.svg
    
elif [ "$1" = "all" ]; then
    rm -rf figs npz logs *.png *.nc *.svg *.pdf

else
    echo "Tell me something!"

fi
