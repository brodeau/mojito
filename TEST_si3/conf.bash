#!/bin/bash

YEAR="1997"

SITRCK_DIR="${HOME}/DEV/sitrack"
MOJITO_DIR="${HOME}/DEV/mojito"

NEMO_CONF="NANUK4"
LIST_NEMO_EXP="BBM00 EVP00"

DATE1="${YEAR}0101"
DATE2="${YEAR}0331"

export SI3DATE1="${YEAR}0101"
export SI3DATE2="${YEAR}0331"

XTRA_SFX_SI3=""

NJPAR=4 ; # number of jobs we can launch in //

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        NJPAR=5
        export DATA_DIR="/MEDIA/data"
        XTRA_SFX_SI3="_LIGHT480"
        LIST_NEMO_EXP="BBM00"
        export iHSS=1 ;  LCOARSEN="640"
        XTRA_SFX_SI3="_LIGHT480"
        DATE2="${YEAR}0115"
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        #
        NJPAR=4       
        #
        ;;
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=30
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac


export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mesh_mask_${NEMO_CONF}_L31_4.2_1stLev.nc"
export FFSM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mask_RGPS_${NEMO_CONF}.nc"

mkdir -p ./figs ./npz

