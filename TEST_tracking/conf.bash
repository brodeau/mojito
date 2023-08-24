#!/bin/bash

YEAR="1997"

SITRCK_DIR="${HOME}/DEV/sitrack"

. ../TEST_rgps/conf.bash

export NEMO_CONF="NANUK4"
LIST_NEMO_EXP="BBM2304 EVP2304"

export SI3DATE1="19961215"
export SI3DATE2="19970420"

XTRA_SFX_SI3=""

ISEED_BASE='selection' ; # 'selection' or 'quads' or 'defs' => will seed based on one or the other

#export DIRIN_PREPARED_RGPS="/home/laurent/tmp/MOJITO/TEST_rgps"
export DIRIN_PREPARED_RGPS="${MOJITO_DIR}/TEST_rgps"


if [ "${MODE}" = "rgps" ]; then
    MODE="rgps_track"
fi

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
        XTRA_SFX_SI3="_LIGHT480"
        LIST_NEMO_EXP="BBM00"
        ISEED_BASE='defs'; # use remaining Quads of deformation...
        #
        ;;
    "frazilo")
        export DATA_DIR="/data"
        #
        #LIST_NEMO_EXP="BBM23A00" ; ISEED_BASE='defs' ; # For scaling
        #LIST_NEMO_EXP="EVP23A00" ; ISEED_BASE='defs' ; # For scaling
        #LIST_NEMO_EXP="BBM23A01" ; ISEED_BASE='defs' ; # For scaling
        LIST_NEMO_EXP="BBM23A02" ; ISEED_BASE='defs' ; # For scaling        
        #
        ;;
    "jackzilla")
        export DATA_DIR="/data1/nobackup/laurent"
        #
        LIST_NEMO_EXP="EVP2305" ; ISEED_BASE='defs' ; # For scaling
        #LIST_NEMO_EXP="BBM2305" ; ISEED_BASE='defs' ; # For scaling
        #LIST_NEMO_EXP="BBM2304" ; ISEED_BASE='defs' ; # For scaling
        #LIST_NEMO_EXP="BBM2304 EVP2304" ; ISEED_BASE='defs' ; # For maps
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

echo; echo " * MIND: mode = ${MODE} !"; echo



export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mesh_mask_${NEMO_CONF}_L31_4.2_1stLev.nc"

mkdir -p ./figs ./npz
