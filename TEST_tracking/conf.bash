#!/bin/bash

YEAR="1997"

MOJITO_DIR="${HOME}/DEV/mojito"

NEMO_CONF="NANUK4"
NEMO_EXP="BBM00"

DT_BINS_H=72  ; # width of a bin for time interpolation (hours)

DATE1="${YEAR}0101"
DATE2="${YEAR}0430"

#LIST_RES="10 20 30"
#LIST_RES="10"
#LIST_STREAM="000 001 002 003 004 005 006 007 008"
#NSTREAMA="0"
#NSTREAMB="8"
MARKER_SIZE=10


FSI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0331_icemod.nc4"

NJPAR=4 ; # number of jobs we can launch in //

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        export DATA_DIR="/MEDIA/data"
        #
        DATE2="${YEAR}0201"
        FNCSEED="RGPS_ice_drift_${YEAR}-01-01_${YEAR}-03-01_lb_LIGHT.nc"
        FSI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0331_icemod_LIGHT480.nc4"
        #
        MARKER_SIZE=20
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        #
        NJPAR=4
        #
        DATE2="${YEAR}0301"
        FNCSEED="RGPS_ice_drift_${YEAR}-01-01_${YEAR}-03-01_lb_LIGHT.nc"
        FSI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0331_icemod_LIGHT480.nc4"
        #
        MARKER_SIZE=20
        ;;
    "frazilo")
        export DATA_DIR="/data"
        #
        NJPAR=28
        #
        DATE2="${YEAR}0401"
        FNCSEED="RGPS_ice_drift_${YEAR}-01-01_${YEAR}-05-01_lb.nc"
        #
        MARKER_SIZE=5
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac



export FNCSEED="${DATA_DIR}/data/mojito/seeding_from_rgps/${FNCSEED}"

export FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_EXP}/${FSI3IN}"

export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mesh_mask_${NEMO_CONF}_L31_4.2_1stLev.nc"

#="
