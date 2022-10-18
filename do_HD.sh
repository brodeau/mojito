#!/bin/bash

if [ $1 -eq 1 ]; then
    ./01_selection_xy.py ${DATA_DIR}/data/RGPS_Kwok_98/RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4 1997
fi


list_streams="000 001 002 003"

if [ $1 -eq 3 ]; then
    
    for cs in ${list_streams}; do

        listF=`\ls npz/SELECTION_buoys_RGPS_stream${cs}_????????_??h??_????????_??h??.npz`

        for ff in ${listF}; do
            ./03_generate_quad_mesh.py ${ff} &
        done

        wait

    done

fi


if [ $1 -eq 4 ]; then

    for cs in ${list_streams}; do
        listQ=`\ls npz/Q-mesh_stream${cs}_????????_??h??_* | head -2`

        echo; echo "./04_deformations.py ${listQ}"

        ./04_deformations.py ${listQ} &
    done

    wait
fi

