#!/bin/bash

. conf.bash

if [ "${1}" = "" ]; then
    echo "USAGE: ${0} <scale> (<orig>)"
    exit
fi

cscal="${1}"

iorig=0
if [ "${2}" != "" ]; then
    iorig=1
    corig="${2}"
else
    lst_exps=`echo ${EXPS} | sed -e s/','/' '/g`
    corig="${lst_exps,,}"
fi

echo; echo " ORIGS: ${corig}"; echo


for cc in ${corig}; do

    echo; echo " *** ${cc}"

    if [ ! -d ./${cc} ]; then
        echo "ERROR: no directory named ./${cc} !"; exit
    fi

    echo " Orig and scale: ${cc} ${cscal}"
    
    rm -f ./bak/${cc}/*_${cscal}km*

    rsync -avP ${cc}/*_${cscal}km* ./bak/${cc}/

    echo
done
