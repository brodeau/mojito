#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

SIZE_MOSAIC=900

#LIST_RES="10 20 30"
LIST_RES="10"

LIST_STREAM="000 001"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        DATA_DIR="/MEDIA/data"
        FILIN="RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4"
        LIST_STREAM="000"
        ;;
    "mcp-oceannext-01")
        DATA_DIR="/data/gcm_setup"
        FILIN="RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4"
        #LIST_STREAM="000 001"
        ;;
    "frazilo")
        DATA_DIR="/data"
        FILIN="RGPS_1996-11-07_1997-06-01_traj.nc4"
        LIST_STREAM="000 001 002 003 004 005"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"

