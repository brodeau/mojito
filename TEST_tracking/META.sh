#!/bin/bash

if [ "${1}" = "0" ]; then
    ./do0CopyRGPSQuads.sh
    wait
fi


./do1Tracking_mp.sh
wait
./do2GenerateQuads_mp.sh
wait
./do3ComputeDef_mp.sh
wait
./do4gatherDef.sh
wait
./do5mkPDFs.sh
wait

