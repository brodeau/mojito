#!/bin/bash

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        RGPS_DIR="/MEDIA/data/data/RGPS_Kwok_98"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac



EXE="${HOME}/DEV/mojito/rgps/01_selection_xy.py"


${EXE} ${RGPS_DIR}/RGPS_1996-11-07_1997-06-01_traj_LIGHT.nc4 1997

