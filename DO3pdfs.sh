#!/bin/bash



RESKM=20
DT_BINS_H=6

mkdir -p figs


for cf in "shear" "divergence" "convergence" "absDiv"; do

    #fRGPS="./TEST_rgps/npz/PDF_RGPS_19970103-19970313_${cf}.npz"
    fRGPS=`find ./TEST_rgps/npz -name PDF_RGPS_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
    
    fBBM=`find ./TEST_tracking/npz -name PDF_NEMO-SI3_NANUK4_BBM00_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
    fEVP=`find ./TEST_tracking/npz -name PDF_NEMO-SI3_NANUK4_EVP00_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
    
    if [ "${fRGPS}" != "" ]; then
        if [ "${fBBM}" != "" ] & [ "${fEVP}" != "" ]; then
            ./plot_comp_PDFs.py ${fRGPS} ${fBBM} ${fEVP}
        elif [ "${fBBM}" != "" ]; then
            ./plot_comp_PDFs.py ${fRGPS} ${fBBM}
        fi
    else
        echo " No ${fRGPS} !"
    fi

             
done
