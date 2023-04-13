#!/bin/bash

#. ../TEST_rgps/conf.bash ; # Get the resolulion "RESKM" !
#. ./conf.bash

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
