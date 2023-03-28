#!/bin/bash

import_frazilo()
{
    rsync -avP -e "ssh -p ${PORT_FRAZILO} -c aes128-ctr" laurent@localhost:${1} ${2}
}

HERE=`pwd`

for cw in "rgps" "tracking"; do

    mkdir -p ${HERE}/TEST_${cw}/npz
    
    cd ${HERE}/TEST_${cw}/npz

    echo; echo " * Inside `pwd` !"
    
    import_frazilo /home/laurent/DEV/mojito/TEST_${cw}/npz/PDF\*.npz .


done
