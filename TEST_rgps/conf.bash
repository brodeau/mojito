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
MODE='rgps'
DEF_EXPORT=''

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
        MODE='rgps'
        #
        # 10km:
        DT_BINS_H=6 ;  RESKM=10
        #DT_BINS_H=72 ;  RESKM=10
        #
        # 20km:
        #DT_BINS_H=6 ;  RESKM=20
        #
        # 40km:
        #DT_BINS_H=6 ; RESKM=40
        #
        # 80km:
        #DT_BINS_H=6 ;  RESKM=80 ; LIST_RD_SS="72 75 78"
        #
        # 160km: in decreasing order of quality:
        #DT_BINS_H=6 ;  RESKM=160 ; LIST_RD_SS="145 150 155 160 165" ; #ok!
        #
        # 320km:        
        #DT_BINS_H=$((24*3)) ;  RESKM=320; LIST_RD_SS="280 285 290 295 300 305 310 313  315  317 320 325 330 335 340 345 350"
        #
        # 640km:
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="570 580 590 600 610 615 620 625 630 633 635 637 640 645 650 655 660     670 680 690 700" ; #WINNA!
        ###DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="570 580 590 600 610 620 625 630 633 635 637 640 645 650 660 670 680 690 700" ; #better than above!                
        #
        ;;        
    "merlat")
        NJPAR=5
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #DT_BINS_H=6 ;         DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10
        DT_BINS_H=6 ;         DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10 ; MODE='rgps' ; DEF_EXPORT='E'
        #DT_BINS_H=$((3*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10 ; MODE='xlose' ; # FOR P90!
        #DT_BINS_H=$((5*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0123"; RESKM=10 ; MODE='xlose' ; # FOR MAPS!!!
        #
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        #FILIN="RGPS_${YEAR}LIGHT.nc4"
        #DT_BINS_H=6 ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10
        DT_BINS_H=$((5*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331"; RESKM=10 ; MODE='xlose' ; # FOR MAPS!!!
        #DT_BINS_H=$((7*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331"; RESKM=10
        #
        #
        ;;
    "jackzilla")
        export DATA_DIR="/data1/nobackup/laurent"
        #
        NJPAR=34
        #
        #DT_BINS_H=6 ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10 ; MODE='xlose' ; # FOR MAPS!!!
        #DT_BINS_H=24 ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10 ; MODE='xlose' ; # FOR MAPS!!!
        #DT_BINS_H=$((3*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0131"; RESKM=10 ; MODE='xlose' ; # FOR MAPS!!!
        DT_BINS_H=$((5*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0331"; RESKM=10 ; MODE='xlose' ; # FOR MAPS!!!
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"


