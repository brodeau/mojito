#!/bin/bash

EXE1="python3 -u /home/laurent/DEV/mojito/util/convert_Qnpz2SeedNC.py"
EXE2="python3 -u /home/laurent/DEV/sitrack/si3_part_tracker.py"
EXE3="python3 -u /home/laurent/DEV/mojito/generate_quad_rccl.py"

rm -f *.png

STR="S004_dt72_19970115"

QUADFILE_MODEL="../TEST_rgps/npz/Q-mesh_RGPS_${STR}-08h44t0_19970115-08h44_10km.npz"

if [ "${1}" = "1" ]; then
    
    ${EXE1} ../TEST_rgps/npz/Q-mesh_RGPS_${STR}
    
    
    fout=`\ls ./nc/Q-mesh_RGPS_${STR}h??_*km.nc`

    
    ${EXE2} -i /MEDIA/data/NANUK4/BBM00/NANUK4_ICE-BBM00_1h_19970101_19970331_icemod_LIGHT480.nc4 \
            -m /MEDIA/data/NANUK4/NANUK4.L31-I/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
            -s  ${fout}
fi



${EXE3} nc/NEMO-SI3_NANUK4_BBM00_tracking_S004_dt72_19970115h10_19970118h20_10km.nc 0,72 \
        ${QUADFILE_MODEL} \
        10 xlose
