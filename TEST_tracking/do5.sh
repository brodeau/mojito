#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/plot_comp_PDFs.py"


for fv in "divergence" "shear"; do

    lst1=`\ls ./npz/PDF_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_${YEAR}????-${YEAR}????_${fv}.npz`
    nbf=`echo ${lst1} | wc -w`
    if [ ${nbf} -ne 1 ]; then echo "PROBLEM: more than 1 file in SI3 stuff!!!"; exit; fi

    lst2=`\ls ../TEST_brgps/npz/PDF_RGPS_S???-S???_${fv}.npz`
    nbf=`echo ${lst2} | wc -w`
    if [ ${nbf} -ne 1 ]; then echo "PROBLEM: more than 1 file in RGPS file!!!"; exit; fi
    
    
    CMD="${EXE} ${lst1} ${lst2}"    
    echo
    echo " *** About to launch:"
    echo "     ${CMD}"
    echo

    ${CMD}

    
done

exit

