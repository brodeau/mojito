#!/bin/bash

if [ "${2}" = "" ]; then
    echo "USAGE: ${0} <orig> <scale>"
    exit
fi

corig="${1}"
cscal="${2}"


if [ ! -d ./${corig} ]; then
    echo "ERROR: no directory named ./${corig} !"; exit
fi
if [ ! -d ./bak/${corig} ]; then
    echo "ERROR: no directory named ./bak/${corig} !"; exit
fi

echo " Orig and scale: ${corig} ${cscal}"




rsync -avP ./bak/${corig}/*_${cscal}km* ${corig}/
