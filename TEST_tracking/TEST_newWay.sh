#!/bin/bash

EXE="python3 -u /home/laurent/DEV/sitrack/si3_part_tracker.py"

${EXE} -i /MEDIA/data/NANUK4/BBM00/NANUK4_ICE-BBM00_1h_19970101_19970331_icemod_LIGHT480.nc4 \
       -m /MEDIA/data/NANUK4/NANUK4.L31-I/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
       -s  ~/DEV/mojito/TEST_rgps/nc/Q-mesh_RGPS_S004_dt72_19970115-08h44t0_19970115-08h44_10km.nc
