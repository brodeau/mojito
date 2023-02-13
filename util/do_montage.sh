#!/bin/bash

# 1131x1080
nx=1131
ny=1080

montage_opt="-background black -geometry $((nx+4))x${ny}+2+0"

list_track=`\ls  track/Pos_buoys_*.png`

mkdir -p ./montage

for ft in ${list_track}; do

    fb=`basename ${ft}`
    date=`echo ${fb} | cut -d_ -f4`

    fr=`\ls rgps/*_${date}*.png`

    fo="./montage/OUT_${date}.png"
    #echo $date $fr

    CMD="montage ${montage_opt} ${fr} ${ft} ${fo}"
    echo
    echo $CMD
    ${CMD}

done


