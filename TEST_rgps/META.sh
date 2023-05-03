#!/bin/bash

. ./conf.bash


echo " *** ${RESKM}km !"


./do1bCoarsify.sh
wait

./do2GenerateQuads_mp.sh
wait
./do3ComputeDef_mp.sh
wait
./do4gatherDef.sh
wait
./do5mkPDFs.sh
wait
