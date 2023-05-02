#!/bin/bash

. ./conf.bash


echo " *** ${RESKM}km !"


#if [ ${RESKM} -gt 10 ]; then
./do1bCoarsify.sh
wait
#fi

./do2GenerateQuads_mp.sh
wait
./do3ComputeDef_mp.sh
wait
./do4gatherDef.sh
wait
./do5mkPDFs.sh
wait
