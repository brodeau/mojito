#!/bin/bash


for dd in 10 15 20; do

    CPREF="npz/Q-mesh_stream000_19970104_18h00"
    
    ../deformation.py ${CPREF}_19970104_18h00_${dd}km_Sampled.npz ${CPREF}_19970107_18h00_${dd}km_Sampled.npz 100

done

