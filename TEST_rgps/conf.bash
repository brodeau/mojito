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
export USE_S='' ; # tells `gatherDeformation` that we need to use the 'SLCT' deformation files...

export NEMO_CONF="NANUK4"

host=`hostname | cut -d '.' -f2`
case ${host} in
    #
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=30
        #
        DATE1="19961215" ; DATE2="${YEAR}0420"
        #
        MODE='rgps' ; DEF_EXPORT='E'
        DT_BINS_H=72
        #
        # 10km:
        #RESKM=10; LIST_RD_SS="5 6 7 8" ; USE_S="S"
        ###RESKM=10 ; USE_S=""
        #
        # 20km:
        #RESKM=20 ; LIST_RD_SS="13 15 17" ; USE_S="S"
        #
        # 40km:
        #RESKM=40 ; LIST_RD_SS="33 35 37 39" ; USE_S="S"
        #
        # 80km:
        #RESKM=80 ; LIST_RD_SS="66 69 72  75  78 81 84 87" ; USE_S="S"
        #
        # 160km:
        RESKM=160 ; LIST_RD_SS="135 140 144 148 152 154  155  156 158 160" ; USE_S="S"
        #
        # 320km:
        ##DT_BINS_H=72 ;  RESKM=320; LIST_RD_SS="290 295 300 305 310 313  315  317 320 325 330 335 340"; export USE_S="S"
        #DT_BINS_H=72 ;  RESKM=320; LIST_RD_SS="270 280 290 295 300 305 310 313   315   317 320 325 330 335 340 350 360"; export USE_S="S"
        #
        # 640km:
        #DT_BINS_H=72 ; RESKM=640 ; LIST_RD_SS="520 530 540 550 560 570 580 590 600 610 615 620 625 630 635 640 645 650 655 660 670 680 690 700 710 720 730 740 750 760"; export USE_S="S"
        #
        ;;        
    "merlat")
        NJPAR=4
        export DATA_DIR="/MEDIA/data"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        #
        #DT_BINS_H=6 ;         DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10 ; MODE='rgps' ; DEF_EXPORT='E'
        #
        #DT_BINS_H=6 ;         DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10        
        #DT_BINS_H=$((3*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10 ; MODE='xlose' ; # FOR P90!
        #DT_BINS_H=$((5*24)) ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0123"; RESKM=10 ; MODE='xlose' ; # FOR MAPS!!!
        #
        DT_BINS_H=6 ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=40 ; MODE='rgps' ; DEF_EXPORT='E' ; LIST_RD_SS="36 38"
        ;;
    "mcp-oceannext-01")
        NJPAR=4
        export DATA_DIR="/data/gcm_setup"
        FILIN="RGPS_${YEAR}LIGHT.nc4"
        DT_BINS_H=6 ; DATE1="${YEAR}0101" ; DATE2="${YEAR}0115"; RESKM=10 ; MODE='rgps' ; DEF_EXPORT='E'
        #
        ;;
    "jackzilla")
        export DATA_DIR="/data1/nobackup/laurent"
        #
        NJPAR=34
        #
        DATE1="19961215" ; DATE2="${YEAR}0420"
        #
        MODE='rgps' ; DEF_EXPORT='E'
        DT_BINS_H=72
        #
        # 10km:
        #RESKM=10; LIST_RD_SS="5 6 7 8" ; USE_S="S"
        #
        # 20km:
        #RESKM=20 ; LIST_RD_SS="13 14 15 16 17" ; USE_S="S"
        #
        # 40km:
        RESKM=40 ; LIST_RD_SS="31 33 35 37 39 41" ; USE_S="S"
        #
        # 80km:
        #RESKM=80 ; LIST_RD_SS="66 69 72 75 78 81 84" ; USE_S="S"; #ok
        #RESKM=80 ; LIST_RD_SS="66 69 72 75 78 81 84 87 90" ; USE_S="S"; #ok
        #
        # 160km: in decreasing order of quality:
        #RESKM=160 ; LIST_RD_SS="144 148 152 156 160 164" ; USE_S="S"
        #
        # 320km:
        ##DT_BINS_H=72 ;  RESKM=320; LIST_RD_SS="260 270 280 285 290 295 300 305 310   315   320 325 330 335 340 345 350 355 360 370 380"; export USE_S="S"
        ####DT_BINS_H=72 ;  RESKM=320; LIST_RD_SS="270 275 280 285 290 295 300 305 310   315   320 325 330 335 340 345 350 355 360"; export USE_S="S"
        #DT_BINS_H=72 ;  RESKM=320; LIST_RD_SS="290 295 300 305 310 313  315  317 320 325 330 335 340"; export USE_S="S"
        # => add more at hihger scales????
        #
        # 640km:
        #DT_BINS_H=72 ; RESKM=640 ; LIST_RD_SS="540 550 560 570 580 590 600 610 615 620 625 630 635 640 645 650 655 660 670 680 690 700 710 720 730 740"; export USE_S="S"
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"


