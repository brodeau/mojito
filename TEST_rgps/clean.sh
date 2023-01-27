#!/bin/bash

if [ "$1" = "2" ]; then
    
    rm -f figs/00_*.png figs/*[TQ]mesh.png npz/[TQ]-mesh_*.npz logs/out_S*__*.out

elif [ "$1" = "3" ]; then

    rm -f figs/z[dst]_* npz/DEFORMATIONS_z* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out

else
    # All!
    rm -rf figs npz *.png
fi
