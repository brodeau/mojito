#!/bin/bash

YEAR="1997"

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
    #
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=30
        #
        DATE1="${YEAR}0101" ; DATE2="${YEAR}0331"
        #
        # 10km:
        #DT_BINS_H=6 ;  RESKM=10
        # 20km:
        #DT_BINS_H=6 ;  RESKM=20 ; # defaut:15
        # 40km:
        #DT_BINS_H=6 ; RESKM=40 ; # defaut:35
        # 80km:
        ##DT_BINS_H=6 ;  RESKM=80 ; LIST_RD_SS="73 83" # defaut:73
        #DT_BINS_H=6 ;  RESKM=80 ;  # defaut:73
        #
        # 160km:
        #DT_BINS_H=6 ;  RESKM=160;  # defaut:145
        ##DT_BINS_H=6 ;  RESKM=160; LIST_RD_SS="130 145 160" # defaut:145
        ##DT_BINS_H=$((24*3)) ;  RESKM=160;  # defaut:145  ; # Ok for at least RGPS!
        ##DT_BINS_H=$((24*3)) ;  RESKM=160; LIST_RD_SS="130 145 160"  ; # defaut:145  ; # Ok for at least RGPS!
        #
        # 320km:
        DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 280 310"  # => p (240, 400) Good!
        #DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="220 310"  # => 77p (240, 400) Good! 
        #DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="310"  # => 47p (240, 400) Good!
        ###DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="310"  #  => 47p (220, 420)
        #
        ##DT_BINS_H=6 ;  RESKM=320; ;  # defaut:295 => 59p
        #
        #DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="310"  # defaut:295 => 47p
        #DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="280 295 310 320 330"  # defaut:295
        ##DT_BINS_H=$((24*3)) ;  RESKM=320;  # defaut:295
        ##DT_BINS_H=$((24*3)) ;  RESKM=320; LIST_RD_SS="280 295 310 325" ; # defaut:295
        #
        # 640km:
        #DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="580 590 595 600 605 610 615 620 625 630 635 640 650" ; # default:620 => 27 points! ok...
        ##DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="580 590 600 610 620 630 640 650 660" ; # default:620 => 25 points! ok...
        ##DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="580 590 600 605 610 615 620 625 630 640 650" ; # default:620 => 27 points! ok...
        

        #
        #
        ;;        
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
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"


