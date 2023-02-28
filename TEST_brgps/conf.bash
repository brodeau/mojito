#!/bin/bash

YEAR="1997"
#YEAR="2007"

MOJITO_DIR="${HOME}/DEV/mojito"

DT_BINS_H=12  ; # width of a bin for time sampling [hours] 

DATE1="${YEAR}0101"
DATE2="${YEAR}0430"

#LIST_RES="10 20 30"
LIST_RES="10"

FILIN="RGPS_${YEAR}.nc4"

MARKER_SIZE=10

NJPAR=4 ; # number of jobs we can launch in //

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        NJPAR=7
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #FILIN="RGPS_1996-12-21_1997-02-01.nc4"
        #
        #DT_BINS_H=6  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=84  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=$((24*7))  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=$((24*10))  ; # width of a bin for time sampling [hours]
        #
        #DATE1="$((YEAR-1))1227_00:00" ; DT_BINS_H=$((24*10))
        DATE1="${YEAR}0103_00:00" ; DT_BINS_H=$((24*10))
        #
        DATE2="${YEAR}0131"
        MARKER_SIZE=10
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        #FILIN="RGPS_${YEAR}LIGHT.nc4"
        FILIN="RGPS_1996-12-21_1997-02-01.nc4"
        #
        #DT_BINS_H=72  ; # width of a bin for time sampling [hours]
        DT_BINS_H=6  ; # width of a bin for time sampling [hours]
        #
        DATE1="${YEAR}0103_00:00" ; DT_BINS_H=$((24*10))
        DATE2="${YEAR}0131"
        #
        MARKER_SIZE=10
        ;;
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=30
        #
        DT_BINS_H=6  ; # width of a bin for time sampling [hours]
        #
        DATE2="${YEAR}0315"
        #
        MARKER_SIZE=50
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"
