#!/bin/bash

if [ "${1}" = "s" ]; then
    ./do0bunchSeed.sh
    ./do1Tracking_mp.sh
fi

if [ "${1}" = "t" ]; then
    ./do1Tracking_mp.sh
fi

./do2GenerateQuads_mp.sh
./do3ComputeDef_mp.sh
./do4gatherDef.sh
./do5mkPDFs.sh


