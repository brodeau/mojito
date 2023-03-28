#!/bin/bash

EXE="./plot_comp_PDFs.py"

RESKM=40
DT_BINS_H=6

mkdir -p figs


for cf in "shear" "divergence" "convergence" "absDiv"; do

    fRGPS=`find ./TEST_rgps/npz -name PDF_RGPS_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
    
    fBBM=`find ./TEST_tracking/npz -name PDF_NEMO-SI3_NANUK4_BBM00_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
    fEVP=`find ./TEST_tracking/npz -name PDF_NEMO-SI3_NANUK4_EVP00_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
    
    if [ "${fRGPS}" != "" ]; then
        if [ "${fBBM}" != "" ] & [ "${fEVP}" != "" ]; then
            ${EXE} ${fRGPS} ${fBBM} ${fEVP}
        elif [ "${fBBM}" != "" ]; then
            ${EXE} ${fRGPS} ${fBBM}
        fi
    else
        echo " No ${fRGPS} !"
    fi

             
done
