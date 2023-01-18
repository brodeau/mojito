#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

DT_BINS_H=12  ; # width of a bin for time sampling [hours] 

SIZE_MOSAIC=900

DATE1="${YEAR}0101"
DATE2="${YEAR}0531"

#LIST_RES="10 20 30"
LIST_RES="10"

FILIN="RGPS_${YEAR}.nc4"

LIST_STREAM="000 001"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        LIST_STREAM="000"
        DATE2="${YEAR}0131"
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        LIST_STREAM="000 001"
        DATE2="${YEAR}0330"
        ;;
    "frazilo")
        export DATA_DIR="/data"
        DT_BINS_H=$((24*7))  ; # width of a bin for time sampling [hours] 
        #LIST_STREAM="000 001 002 003 004 005"
        LIST_STREAM="000 001"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"

