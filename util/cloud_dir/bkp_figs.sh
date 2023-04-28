#!/bin/bash


#if [ "${1}" = "" ]; then
#    echo "USAGE: ${0} <scale>"
#    exit
#fi


listo=`find ./figs -name "*_dt0*_o.png" | grep -v '_oo.png'`
list=`find  ./figs -name "*_dt0.png"    | grep -v   'o.png'`

echo $list ; echo; sleep 1

for ff in ${listo}; do

    fn=`echo ${ff} | sed -e "s|dt0_o.png|dt0_oo.png|g"`

    
    echo "mv -f ${ff} ${fn}"
    
    mv -f ${ff} ${fn}

done


for ff in ${list}; do

    fn=`echo ${ff} | sed -e "s|dt0.png|dt0_o.png|g"`

    
    echo "mv -f ${ff} ${fn}"
    
    mv -f ${ff} ${fn}

done
