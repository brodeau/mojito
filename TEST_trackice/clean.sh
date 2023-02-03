#!/bin/bash

rm -f mesh_mask*.nc *.svg

. conf.bash

if [ "$1" = "1" ]; then
    
    rm -f npz/[TQ]-mesh_*.npz logs/out_[QT]-mesh*_*.out logs/err_[QT]-mesh*_*.err

elif [ "$1" = "2" ]; then
    rm -f figs/z[dst]_* npz/DEFORMATIONS_${NEMOCONF}_*
    rm -f logs/${NEMOCONF}_${EXPRMNT}*.out logs/${NEMOCONF}_${EXPRMNT}*.err
    
elif [ "$1" = "all" ]; then
    rm -rf figs npz logs *.png

else
    echo "Tell me something!"

fi
