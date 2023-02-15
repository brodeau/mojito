#!/bin/bash

if [ "$1" = "2" ]; then
    
    rm -f figs/00_*.png figs/*[TQ]mesh.png npz/[TQ]-mesh_*.npz logs/out_S*__*.out figs/fig0*.png

elif [ "$1" = "3" ]; then
    killall -9 deformation.py do3.sh
    rm -f figs/deformation/*.png npz/DEFORMATIONS_* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out
    
elif [ "$1" = "4" ]; then
    rm -f npz/PDF_*.npz *PDF*.svg
    
elif [ "$1" = "all" ]; then
    rm -rf figs npz logs *.png *.nc *.svg *.pdf

else
    echo "Tell me something!"

fi
