#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"
DT_BINS_H=12  ; # width of a bin for time sampling [hours] 
DATE1="${YEAR}0101"
DATE2="${YEAR}0430"
FILIN="RGPS_${YEAR}.nc4"
NJPAR=4 ; # number of jobs we can launch in //
LIST_RD_SS=""
LIST_MINDC=""
RESKM=10

host=`hostname | cut -d '.' -f2`
case ${host} in
    #
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=30
        #
        DATE1="${YEAR}0101" ; DATE2="${YEAR}0331"
        #
        # 10km:
        #DT_BINS_H=$((24*3)); RESKM=10
        #DT_BINS_H=6 ;  RESKM=10
        #
        # 20km:
        #DT_BINS_H=$((24*3)); RESKM=20
        #
        # 40km:
        #DT_BINS_H=$((24*3)); RESKM=40
        #
        # 80km:
        #DT_BINS_H=$((24*3)); RESKM=80
        #
        # 160km:
        #DT_BINS_H=$((24*3)); RESKM=160
        #
        # 320km:
        #DT_BINS_H=$((24*3)); RESKM=320; LIST_RD_SS="270 295 320"
        #
        # 640km:
        DT_BINS_H=$((24*3)); RESKM=640; LIST_RD_SS="590 620 650"
        #
        ;;        
    "merlat")
        NJPAR=5
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        DT_BINS_H=6 ;         DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=10
        #DT_BINS_H=6 ;         DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=40
        #DT_BINS_H=$((3*24)) ; DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0131"; RESKM=160
        #
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #FILIN="RGPS_1996-12-21_1997-02-01.nc4"
        #
        #DT_BINS_H=72  ; # width of a bin for time sampling [hours]
        #DT_BINS_H=6  ; DATE2="${YEAR}0131"
        #
        #DATE1="${YEAR}0103_00:00" ; DT_BINS_H=$((24*5)) ; DATE2="${YEAR}0115"
        #
        FILIN="RGPS_${YEAR}.nc4"; DATE1="${YEAR}0103_00:00" ; DT_BINS_H=$((24*7)) ; DATE2="${YEAR}0109" ; # Best display of full coverage ?
        
        #DATE1="${YEAR}0101_00:00" ; DT_BINS_H=$((24*9)) ; DATE2="${YEAR}0120" ; # testing for 200 km
        #
        #DATE1="${YEAR}0103_00:00" ; DATE2="${YEAR}0119"; DT_BINS_H=$((24*7))
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"


