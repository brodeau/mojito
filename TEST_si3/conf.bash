#!/bin/bash

#YEAR="1997"
YEAR="2016"

SITRCK_DIR="${HOME}/DEV/sitrack"
MOJITO_DIR="${HOME}/DEV/mojito"

export NEMO_CONF="HUDSON4"; NLEV=16
LIST_NEMO_EXP="BBM00"

DATE1="${YEAR}0101"
DATE2="${YEAR}0331"

export LCOARSEN=( "10" "20" "40" "80" "160" "320" "640" )
export LDTINCRM=(  "3"  "3"  "3"  "3"  "3"   "3"   "1"  )

export SI3DATE1="${YEAR}0101"
export SI3DATE2="${YEAR}0331"


export FMASK_RGPS="tmask_HUDSON4.nc"

MODE="model"

XTRA_SFX_SI3=""

NJPAR=4 ; # number of jobs we can launch in //

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        NJPAR=5
        export DATA_DIR="/MEDIA/data"
        XTRA_SFX_SI3="_LIGHT480"
        LIST_NEMO_EXP="BBM00"
        export iHSS=1 ;  LCOARSEN="320"
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
        export LCOARSEN=( "10" ); export LDTINCRM=( "3" )
        export MODE="rgps_map" ; #xlose
        #
        #LIST_NEMO_EXP="BBM2302 EVP2302"
        #LIST_NEMO_EXP="EVP2302"
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L${NLEV}-I/mesh_mask_${NEMO_CONF}_L${NLEV}_4.2_1stLev.nc"
export FFSM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L${NLEV}-I/${FMASK_RGPS}"

mkdir -p ./figs ./npz

