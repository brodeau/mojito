#!/bin/bash

YEAR="1997"

SITRCK_DIR="${HOME}/DEV/sitrack"

. ../TEST_rgps/conf.bash



NEMO_CONF="NANUK4"
LIST_NEMO_EXP="BBM2302 EVP2302"

export SI3DATE1="${YEAR}0101"
export SI3DATE2="${YEAR}0331"

XTRA_SFX_SI3=""

ISEED_BASE='selection' ; # 'selection' or 'quads' or 'defs' => will seed based on one or the other

#export DIRIN_PREPARED_RGPS="/home/laurent/tmp/MOJITO/TEST_rgps"
export DIRIN_PREPARED_RGPS="${MOJITO_DIR}/TEST_rgps"


host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        export DATA_DIR="/MEDIA/data"
        XTRA_SFX_SI3="_LIGHT480"
        LIST_NEMO_EXP="BBM00"
        #ISEED_BASE='quads'
        ISEED_BASE='defs'; # use remaining Quads of deformation...
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        LIST_NEMO_EXP="BBM00"
        ISEED_BASE='defs'; # use remaining Quads of deformation...
        #
        ;;
    "frazilo")
        export DATA_DIR="/data"
        #
        ISEED_BASE='defs'
        #
        # For maps:
        #LIST_NEMO_EXP="BBM00"
        #LIST_NEMO_EXP="BBM2302"
        #LIST_NEMO_EXP="EVP2302"
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

echo " * MIND: mode = ${MODE} !"
sleep 1


export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mesh_mask_${NEMO_CONF}_L31_4.2_1stLev.nc"

mkdir -p ./figs ./npz

