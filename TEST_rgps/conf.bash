#!/bin/bash

YEAR="1997"
#YEAR="2007"

MOJITO_DIR="${HOME}/DEV/mojito"

DT_BINS_H=12  ; # width of a bin for time sampling [hours] 

DATE1="${YEAR}0101"
DATE2="${YEAR}0430"

FILIN="RGPS_${YEAR}.nc4"

NJPAR=4 ; # number of jobs we can launch in //


RESKM=10

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
        #DATE1="${YEAR}0103_00:00" ; DT_BINS_H=$((24*10))
        #DATE1="${YEAR}0103_00:00" ; DT_BINS_H=$((24*5)) ; DATE2="${YEAR}0115"
        DATE1="${YEAR}0103_00:00" ; DT_BINS_H=24 ; DATE2="${YEAR}0115" ; RESKM=10
        #DT_BINS_H=6 ; DATE2="${YEAR}0131"
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
        # 10km:
        #DT_BINS_H=6 ; DATE2="${YEAR}0315" ; RESKM=10 ; #scale 10km
        DT_BINS_H=$((24*5)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=10
        # 20km:
        #DT_BINS_H=6 ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=20 ; # => 47806  points for shear !!!
        # 40km:
        #DT_BINS_H=6 ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=40 ; # =>   points for shear !!!
        # 80km:
        #DT_BINS_H=$((24*3)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=80 ; # => 3913 points for shear !!!
        # 160km:
        #DT_BINS_H=$((24*3)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=160 ; # =>  329   points for shear !!!        
        #DT_BINS_H=$((24*5)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=160 ; # =>  Less!!! points for shear !!!
        # 320km:
        #DT_BINS_H=$((24*5)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=320 ; # =>    points for shear !!!
        # 640km:
        #DT_BINS_H=$((24*5)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331" ; RESKM=640 ; # =>    points for shear !!!
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"
