#!/bin/bash

. ./conf.bash

echo ; echo " *** ${RESKM}km !"; echo

if [ "${1}" = "0" ]; then
    ./do0ScanBatches.sh
    wait
fi

./do1Coarsify.sh
wait

./do2GenerateQuads_mp.sh
wait
./do3ComputeDef_mp.sh
wait
./do4gatherDef.sh
wait
./do5mkPDFs.sh
wait
