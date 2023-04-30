#!/bin/bash

if [ "${1}" = "s" ]; then
    ./do0bunchSeed.sh
fi

./do1Tracking_mp.sh
./do2GenerateQuads_mp.sh
./do3ComputeDef_mp.sh
./do4gatherDef.sh
./do5mkPDFs.sh


