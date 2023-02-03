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

#LIST_STREAM="000 001 002 003 004 005 006 007 008"
NSTREAMA="0"
NSTREAMB="8"

MARKER_SIZE=10

NJPAR=4 ; # number of jobs we can launch in //

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #FILIN="RGPS_${YEAR}.nc4"
        DT_BINS_H=6  ; # width of a bin for time sampling [hours]
        #DATE1="$((YEAR-1))1230"        
        DATE2="${YEAR}0201"
        #LIST_RES="20"
        MARKER_SIZE=40
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #LIST_STREAM="000 001 002 003 004 005 006 008 009 010 011"
        NSTREAMB="11"
        #
        DATE1="${YEAR}0101"
        DATE2="${YEAR}0201"
        DT_BINS_H=6  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=$((3*24))  ; # width of a bin for time sampling [hours]
        ;;
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=28
        #
        #DT_BINS_H=$((24*7))  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=$((12))  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=48  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=24  ; # width of a bin for time sampling [hours]
        DT_BINS_H=6  ; # width of a bin for time sampling [hours]
        #
        #NSTREAMB="1"        
        NSTREAMB="136"
        #DATE1="$((YEAR-1))1230"
        #DATE2="${YEAR}0315"
        #
        DATE1="${YEAR}0101"
        #DATE2="${YEAR}0115"
        DATE2="${YEAR}0401"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"

