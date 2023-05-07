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
        #DT_BINS_H=6 ;  RESKM=10
        #
        # 20km:
        #DT_BINS_H=6 ;  RESKM=20
        #
        # 40km:
        #DT_BINS_H=6 ; RESKM=40
        #
        # 80km:
        #DT_BINS_H=6 ;  RESKM=80
        #
        # 160km: in decreasing order of quality:
        #DT_BINS_H=$((24*3)) ;  RESKM=160 ; LIST_RD_SS="145 155 165"
        ##DT_BINS_H=$((24*3)) ;  RESKM=160 ; # 541p        
        ##DT_BINS_H=6 ;  RESKM=160 ; # 571p
        #
        #
        # 320km:
        #DT_BINS_H=$((24*3)) ;  RESKM=320; LIST_RD_SS="290 300 305 310 315 320 325 330 340" ; # The cream!!!
        #
        #
        # 640km:
        DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="590 600 610 620 625 628 631 633 635 637 639 642 645 650 660 670 680"  ; # 
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="590 600 610 620 625 630 633 635 637 640 645 650 660 670 680"  ; # 
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="615 625 630 635 640 645 655"  ; # Pretty good!!! :)
        ##DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="615 625 635 645 655"  ; # ok !!!
        ##DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="600 620 635 650 670"
        #
        ;;        
    "merlat")
        NJPAR=5
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #DT_BINS_H=6 ;         DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=10
        DT_BINS_H=6 ;         DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=40
        #DT_BINS_H=6 ;         DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=80
        #DT_BINS_H=$((3*24)) ; DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0131"; RESKM=160
        #
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        DT_BINS_H=6 ; DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=10
        #DT_BINS_H=6 ; DATE1="${YEAR}0101_00:00" ; DATE2="${YEAR}0115"; RESKM=20
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"


