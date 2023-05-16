#!/bin/bash


#if [ "${1}" = "" ]; then
#    echo "USAGE: ${0} <scale>"
#    exit
#fi


listo=`find ./figs -name "*_o.png"`
list=`find  ./figs -name "*.png" | grep -v 'o.png'`

echo $list ; echo; sleep 1

for ff in ${listo}; do

    fn=`echo ${ff} | sed -e "s|_o.png|_oo.png|g"`

    
    echo "mv -f ${ff} ${fn}"
    
    mv -f ${ff} ${fn}

done


for ff in ${list}; do

    fn=`echo ${ff} | sed -e "s|.png|_o.png|g"`
    
    echo "mv -f ${ff} ${fn}"
    
    mv -f ${ff} ${fn}

done
