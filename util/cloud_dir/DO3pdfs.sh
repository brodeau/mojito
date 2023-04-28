#!/bin/bash

EXE="${HOME}/DEV/mojito/plot_comp_PDFs.py"


#VRESKM=( 10 20 40 80 160 320 640 )
#VDTBIN=(  6  6  6 72  72  72  72 )
VRESKM=(  320 640 )
VDTBIN=(   72  72 )

mkdir -p figs

ic=0
for RESKM in ${VRESKM[*]}; do

    DT_BINS_H=${VDTBIN[${ic}]}
    #RESKM=10 ; DT_BINS_H=6
    #RESKM=20 ; DT_BINS_H=6
    #RESKM=40 ; DT_BINS_H=6
    #RESKM=80 ; DT_BINS_H=72
    #RESKM=160 ; DT_BINS_H=72
    #RESKM=320 ; DT_BINS_H=120
    #RESKM=640 ; DT_BINS_H=120


    echo; echo " *** RESKM=${RESKM} ; DT_BINS_H=${DT_BINS_H} !"


    for cf in "shear" "total" "divergence" "convergence" "absDiv"; do

        fRGPS=`find ./rgps -name PDF_RGPS_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
        echo " * fRGPS = ${fRGPS}"; echo

        fBBM=`find ./bbm -name PDF_NEMO-SI3_NANUK4_BBM00_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
        echo " * fBBM = ${fBBM}"; echo

        fEVP=`find ./evp -name PDF_NEMO-SI3_NANUK4_EVP00_dt${DT_BINS_H}_${RESKM}km_*_${cf}.npz 2>/dev/null`
        echo " * fEVP = ${fEVP}"; echo

        if [ "${fRGPS}" != "" ]; then
            if [ "${fBBM}" != "" ] & [ "${fEVP}" != "" ]; then
                CMD="${EXE} ${fRGPS} ${fBBM} ${fEVP}"
            elif [ "${fBBM}" != "" ]; then
                CMD="${EXE} ${fRGPS} ${fBBM}"
            fi
        else
            echo " No ${fRGPS} !"
        fi

        echo
        echo " *** Wil launch: ${CMD}"
        echo
        ${CMD}
        echo


    done

    ic=`expr ${ic} + 1`
done
