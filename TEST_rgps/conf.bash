#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        RGPS_DIR="/MEDIA/data/data/RGPS_Kwok_98"
        FILIN="RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4"
        ;;
    "frazilo")
        RGPS_DIR="/data/data/RGPS_Kwok_98"
        FILIN="RGPS_1996-11-07_1997-06-01_traj.nc4"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac


