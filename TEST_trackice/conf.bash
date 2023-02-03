#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

NEMOCONF="NANUK4"
EXPRMNT="ICE-BBM00_1h"

HSS=1 ; RESOL0=12 ; MARKER_SIZE=100 ; # NANUK4
TSS=72

#LIST_RES="10 20 30"
#LIST_RES="10"
#LIST_STREAM="000 001"

NbRec=11

NJPAR=4 ; # number of jobs we can launch in //

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        DATA_DIR="/MEDIA/data"
        HSS=3 ; RESOL0=36 ; MARKER_SIZE=100 ; # NANUK4 with sub-samp of 5 !
        ;;
    "mcp-oceannext-01")
        DATA_DIR="/data/gcm_setup"
        HSS=5 ; RESOL0=60 ; MARKER_SIZE=250 ; # NANUK4 with sub-samp of 5 !
        ;;
    "frazilo")
        DATA_DIR="/data"
        NJPAR=28
        #HSS=1 ; RESOL0=12; NbRec=11; MARKER_SIZE=8 ; # NANUK4 no subsampling
        HSS=2 ; RESOL0=25; NbRec=30; MARKER_SIZE=24 ; # NANUK4 no subsampling
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

MESH_MASK="${DATA_DIR}/NANUK4/NANUK4.L31-I/mesh_mask_NANUK4_L31_4.2.nc"

FILIN="${DATA_DIR}/data/trackice_output/${NEMOCONF}/${NEMOCONF}_${EXPRMNT}_${YEAR}0101_${YEAR}0331_HSS${HSS}_TSS${TSS}.npz"


for ff in ${MESH_MASK} ${FILIN}; do
    if [ ! -f ${ff} ]; then
        echo " ${ff} is missing!!!"
        exit
    fi
done
