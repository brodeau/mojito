#!/bin/bash

. ./conf.bash

fin=`\ls ./nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking_nemoTsi3_${DATE1}h00_${YEAR}????h??.nc`


if [ `echo ${fin} | wc -w ` -ne 1 ]; then
    echo "ERROR: problem with available NC files!!!"
    exit
fi





../../generate_quad_mesh.py ${fin}  0,72   ${RESKM}

../../generate_quad_mesh.py ${fin}  72,144 ${RESKM}
