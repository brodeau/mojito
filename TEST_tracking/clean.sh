#!/bin/bash

#if [ "$1" = "2" ]; then
    
#    rm -f figs/00_*.png figs/*[TQ]mesh.png npz/[TQ]-mesh_*.npz logs/out_S*__*.out

#elif [ "$1" = "3" ]; then
#    killall -9 deformation.py do3.sh
#    rm -f figs/z*.png npz/DEFORMATIONS_z* logs/err_Q-mesh_*.err logs/out_Q-mesh_*.out
    
#elif [ "$1" = "all" ]; then
rm -rf figs npz logs *.png *.nc out

#else
#    echo "Tell me something!"

#fi
