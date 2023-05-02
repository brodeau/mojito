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
        DT_BINS_H=6 ;  RESKM=10
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
        # 160km:
        #DT_BINS_H=6 ;  RESKM=160
        #DT_BINS_H=$((24*3)) ;  RESKM=160
        #
        # 320km:
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="300" ; #(260, 380) / 128 pt / Pas mal!
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="295" ; #(260, 380) /  146 pt / Pas mal+!
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="285 295 305" ; #(260, 380) /  444 pt / Pas mal++! ZOLEIL *_*
        #DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="285 295 305" ; #(280, 360) / ZOLEIL++  8==D 
        ##
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="290" ; #(260, 380) / 144 pt / drop!!!
        ##
        ##DT_BINS_H=$((24*3)) ;  RESKM=320; LIST_RD_SS="300" ; #(260, 380) / 124 pt similar
        ##
        #DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 290 300 310 320 335 350"  # =>  266p (260, 380) Good!!!
        #DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 290 295 300 310 320"  # =>  259p (260, 380) pas mal!!!
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 290 295 300 305 310 320 330"  # =>  298p (260, 380) moins bien!
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 290 295 300 310"  # =>  240p (260, 380) ok!
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 290 295 300 310"  # =>  280p (255, 385) ok!
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 290 295 300 310"  # => 236p (260, 380) bof
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 295 310"  # => 167p (260, 380) bof
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 265 280 295 310"  # => 263p (240, 400) Better!
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="250 280 310"  # => 158p (240, 400) Good!
        ##
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="310"  # defaut:295 => 47p
        ##DT_BINS_H=6 ;  RESKM=320; LIST_RD_SS="280 295 310 320 330"  # defaut:295
        #
        # 640km:
        #***********************
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="610 630 650" ; # (480, 800)/Tdev=24h =>  192p / Au dessus c l'soleil! FOLEIL
        #***********************
        #
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="470 500 540 550 560 570 580 590 595 600 605 610 620 630 640 650 660" ; # (480, 800)/Tdev=24h => 1270 points GOOD!!!
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="540 550 560 570 580 590 600 610 620 630 640 650 660" ; # (480, 800)/Tdev=24h => 962 points The Best so far *_*
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="570 580 590 600 610 620 630 640 650 660" ; # (480, 800)/Tdev=24h => 714 p / slight improvement
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="590 600 610 620 630 640 650 660" ; # (480, 800)/Tdev=24h => 549 p / slight improvement
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="580 600 620 640 660" ; # (480, 800)/Tdev=24h =>  352 p / DA BEST!
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="590 620 650" ; # (480, 800)/Tdev=24h =>  211 p / bit less good
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="600 620 640" ; # (480, 800)/Tdev=24h =>  209 p / DA BEST! (mainly due to improvement q1!)
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="610 630 650" ; # (480, 800)/Tdev=24h =>  192p / Au dessus c l'soleil! FOLEIL
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="610 630 650" ; # (520, 760)/Tdev=24h =>   BIN 177p / bit less good
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="615 635 655" ; # (480, 800)/Tdev=24h =>   201p / Au dessus c l'soleil!v3
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="600 620 640" ; # (480, 800)/Tdev=24h =>  209p / Au dessus c l'soleil!v2
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="620 640 660" ; # (480, 800)/Tdev=24h =>   194p / moins bon pur q1 surtout
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="615 630 645" ; # (480, 800)/Tdev=24h =>  205p / bcp moins bien
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="630" ; # (480, 800)/Tdev=24h => pas bon
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="620 640" ; # (480, 800)/Tdev=24h => 
        #
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="470 500 540 550 560 570 580 590 595 600 605 610 620 630 640 650 660" ; # (500, 780)/Tdev=24h =>  1156 points GOOD! (better q1,less for q3)
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="470 500 540 550 560 570 580 590 595 600 605 610 620 630 640 650 660" ; # (520, 760)/Tdev=24h =>  979 points GOOD! (better q1,less for q3)
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="500 540 560 580 600 620 640 660 700" ; # (520, 760)/Tdev=24h =>  571 points good! (better for q2,q3,less q1)
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="500 540 560 580 600 620 640 660 700" ; # (500, 780)/Tdev=24h =>   623 points good!
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="510 540 560 580 600 620 640 660 690 720" ; # (520, 760)/Tdev=24h => 598 points GOOD!
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="520 540 560 580 600 620 640 660 680 700 720 740" ; # (520, 760)/Tdev=24h => 685 points less GOOD!
        #
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="540 560 580 600 610 620 630 640 660 680 700 720 740" ; # (520, 760)/Tdev=24h =>  723 points more GOOD!
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="540 560 580 600 610 620 630 640 660 680 700 720 740" ; # (600, 680)/Tdev=24h =>  346 points trop remonte!
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="540 560 580 600 610 620 630 640 660 680 700 720 740" ; # (500, 780)/Tdev=24h =>  769 points slightly betta!

        
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="500 540 560 580 600 620 640 660 700" ; # (500, 780)/Tdev=18h =>   bad
        #
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="470 500 540 550 560 570 580 590 595 600 605 610 620 630 640 650 660" ; # (500, 780)/Tdev=12h => 254 points bad!
        #DT_BINS_H=$((24*3)) ; RESKM=640 ; LIST_RD_SS="470 500 540 550 560 570 580 590 595 600 605 610 620 630 640 650 660" ; # (480, 800)/Tdev=12h => 296 points bad!
        #
        ##LIST_RD_SS="500 520 530 540 550 560 570 580 590 600 610 620 630 640 650" ; # 134p (520, 760) |
        ##LIST_RD_SS="500 520 540 560 570 575 580 585 590 595 600 605 610 615 620 630" ; # 151p (510, 770) | pas mal
        ##LIST_RD_SS="500 520 540 560 570 580 590 600 610 620 630" ; # Good! 105p (510, 770)
        ##LIST_MINDC="300 200 400 300 200 450 250 400 200 450 250"
        #DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="500 520 540 560 570 580 590 600 610 620 630" ; # default:620 ;  100 points!  (520, 760)
        #DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="500 520 560 570 580 590 600 610 620 630" ; # default:620 ; # The best so far !!!
        #DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="500 520 560 580 590 600 610 620 630" ; # default:620 ; # even better :D
        #DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="520 560 580 590 600 610 620 630" ; # default:620 ; # Better :)
        #DT_BINS_H=$((24*3)) ; RESKM=640; LIST_RD_SS="520 560 580 590 600 610 620" ; # default:620 ; # Not bad!       
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


