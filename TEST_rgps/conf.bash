#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

SIZE_MOSAIC=900

DATE1="${YEAR}0101"
DATE2="${YEAR}0531"

#LIST_RES="10 20 30"
LIST_RES="10"

LIST_STREAM="000 001"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4"
        #FILIN="RGPS_1996-11-07_1997-06-01_traj.nc4"
        LIST_STREAM="000"
        DATE2="${YEAR}0131"
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        FILIN="RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4"
        LIST_STREAM="000 001"
        DATE2="${YEAR}0330"
        ;;
    "frazilo")
        export DATA_DIR="/data"
        FILIN="RGPS_1996-11-07_1997-06-01_traj.nc4"
        #LIST_STREAM="000 001 002 003 004 005"
        LIST_STREAM="000 001"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"

