#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

NEMO_CONF="NANUK4"
NEMO_EXP="BBM00"

DATE1="${YEAR}0101"

FSI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0331_icemod.nc4"

NJPAR=4 ; # number of jobs we can launch in //

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        export DATA_DIR="/MEDIA/data"
        #
        FSI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0331_icemod_LIGHT480.nc4"
        #
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        #
        NJPAR=4       
        #
        #FSI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0331_icemod_LIGHT480.nc4"
        #
        ;;
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=30
        #
        #SI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0205_icemod.nc4" ; # 1 month !!!
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac


DIRIN_PREPARED_RGPS="${MOJITO_DIR}/TEST_rgps/nc"


export FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_EXP}/${FSI3IN}"

export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mesh_mask_${NEMO_CONF}_L31_4.2_1stLev.nc"

mkdir -p ./figs ./npz

