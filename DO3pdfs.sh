#!/bin/bash







for cf in "shear" "divergence" "convergence" "absDiv"; do

    fRGPS="./TEST_rgps/npz/PDF_RGPS_19970103-19970313_${cf}.npz"
    fBBM="./TEST_tracking/npz/PDF_NEMO-SI3_NANUK4_BBM00_19970103-19970313_${cf}.npz"
    fEVP="./TEST_tracking/npz/PDF_NEMO-SI3_NANUK4_EVP00_19970103-19970313_${cf}.npz"

    ./plot_comp_PDFs.py ${fRGPS} ${fBBM} ${fEVP}

done
