#!/bin/bash

YEAR="1997"
#YEAR="2007"

MOJITO_DIR="${HOME}/DEV/mojito"

DT_BINS_H=12  ; # width of a bin for time sampling [hours] 

DATE1="${YEAR}0101"
DATE2="${YEAR}0430"

FILIN="RGPS_${YEAR}.nc4"

NJPAR=4 ; # number of jobs we can launch in //

LIST_RD_SS=""

RESKM=10

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        NJPAR=5
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #DT_BINS_H=6 ;         DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=10
        DT_BINS_H=$((3*24)) ; DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=10
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
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=30
        #
        DATE1="${YEAR}0101" ; DATE2="${YEAR}0331"
        #
        # 10km:
        #DT_BINS_H=6 ;  RESKM=10
        DT_BINS_H=$((24*3)) ;  RESKM=10
        # 20km:
        #DT_BINS_H=6 ;  RESKM=20 ; LIST_RD_SS="15 17 19"
        #DT_BINS_H=$((24*3)) ;  RESKM=20 ; LIST_RD_SS="15 17 19"
        # 40km:
        #DT_BINS_H=6 ;  RESKM=40 ; LIST_RD_SS="30 35 40"
        #DT_BINS_H=$((24*3)) ;  RESKM=40 ; LIST_RD_SS="30 35 40"
        # 80km:
        #DT_BINS_H=$((24*3)) ;  RESKM=80 ; LIST_RD_SS="65 70 75 80 85"
        # 160km:
        #DT_BINS_H=$((24*3)) ;  RESKM=160; LIST_RD_SS="120 130 140 150 160 170 180"
        # 320km:
        #DT_BINS_H=$((24*3)) ;  RESKM=320; LIST_RD_SS="230 240 260 270 280 285 290 295 300 305 310 315 320 325 330 335 340 350"
        # 640km:
        #DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="550 560 570 580 590 600 610 615 620 625 630 635 640 645 650 660 670 680 690"
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"


