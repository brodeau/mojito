#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

#LIST_RES="10 15 20"
LIST_RES="10"

LIST_STREAM="000 001"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        RGPS_DIR="/MEDIA/data/data/RGPS_Kwok_98"
        FILIN="RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4"
        LIST_STREAM="000 001"
        ;;
    "frazilo")
        RGPS_DIR="/data/data/RGPS_Kwok_98"
        FILIN="RGPS_1996-11-07_1997-06-01_traj.nc4"
        LIST_STREAM="001"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac


