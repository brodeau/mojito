#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

SIZE_MOSAIC=900

#LIST_RES="10 20 30"
LIST_RES="10"

LIST_STREAM="000 001"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        DATA_DIR="/MEDIA/data"
        FILIN="./light/NANUK4_ICE-BBM00_1h_19970101_19970331_HSS5_TSS72.npz"
        ;;
    "mcp-oceannext-01")
        DATA_DIR="/data/gcm_setup"
        FILIN="./light/NANUK4_ICE-BBM00_1h_19970101_19970331_HSS5_TSS72.npz"
        ;;
    "frazilo")
        DATA_DIR="/data"
        FILIN="./hd/NANUK4_ICE-BBM00_1h_19970101_19970331_HSS1_TSS72.npz"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

#export RGPS_DIR="${DATA_DIR}/data/RGPS_Kwok_98"



MESH_MASK="${DATA_DIR}/NANUK4/NANUK4.L31-I/mesh_mask_NANUK4_L31_4.2.nc"
