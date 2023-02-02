#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

NEMOCONF="NANUK4"

HSS=1 ; RESOL0=12 ; MARKER_SIZE=100 ; # NANUK4
TSS=72

#LIST_RES="10 20 30"
#LIST_RES="10"
#LIST_STREAM="000 001"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        DATA_DIR="/MEDIA/data"
        HSS=5 ; RESOL0=60 ; MARKER_SIZE=250 ; # NANUK4 with sub-samp of 5 !
        ;;
    "mcp-oceannext-01")
        DATA_DIR="/data/gcm_setup"
        HSS=5 ; RESOL0=60 ; MARKER_SIZE=250 ; # NANUK4 with sub-samp of 5 !
        ;;
    "frazilo")
        DATA_DIR="/data"
        HSS=1 ; RESOL0=12; MARKER_SIZE=5 ; # NANUK4 no subsampling
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

MESH_MASK="${DATA_DIR}/NANUK4/NANUK4.L31-I/mesh_mask_NANUK4_L31_4.2.nc"

FILIN="${DATA_DIR}/data/trackice_output/${NEMOCONF}/${NEMOCONF}_ICE-BBM00_1h_19970101_19970331_HSS${HSS}_TSS${TSS}.npz"


for ff in ${MESH_MASK} ${FILIN}; do
    if [ ! -f ${ff} ]; then
        echo " ${ff} is missing!!!"
        exit
    fi
done
