#!/bin/bash


if [ "${1}" = "" ]; then
    echo "USAGE: ${0} <scale>"
    exit
fi


list=`find . -name "*_${1}km_*.npz" | grep -v bak`

echo $list ; echo; sleep 4

for ff in ${list}; do

    echo "rm ${ff}"
    rm ${ff}

done
