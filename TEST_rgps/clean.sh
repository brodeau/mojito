#!/bin/bash

if [ $1 -eq 2 ]; then
    
    rm -f figs/00_*.png figs/*[TQ]mesh.png npz/[TQ]-mesh_*.npz

else
    # All!
    rm -rf figs npz *.png
fi
