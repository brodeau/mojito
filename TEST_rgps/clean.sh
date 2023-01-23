#!/bin/bash

if [ "$1" = "2" ]; then
    
    rm -f figs/00_*.png figs/*[TQ]mesh.png npz/[TQ]-mesh_*.npz

elif [ "$1" = "3" ]; then

    rm -f figs/z_* npz/DEFORMATIONS_z_*

else
    # All!
    rm -rf figs npz *.png
fi
