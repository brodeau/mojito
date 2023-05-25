#!/bin/bash

EXE1="python3 -u /home/laurent/DEV/mojito/util/convert_Qnpz2SeedNC.py"
EXE2="python3 -u /home/laurent/DEV/sitrack/si3_part_tracker.py"

rm -f *.png

STR="S004_dt72_19970115-08h44t0"

${EXE1} ../TEST_rgps/npz//Q-mesh_RGPS_${STR}


${EXE2} -i /MEDIA/data/NANUK4/BBM00/NANUK4_ICE-BBM00_1h_19970101_19970331_icemod_LIGHT480.nc4 \
       -m /MEDIA/data/NANUK4/NANUK4.L31-I/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
       -s  ~/DEV/mojito/TEST_rgps/nc/Q-mesh_RGPS_${STR}.nc
