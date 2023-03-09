#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/coarsify_point_cloud.py"

# Populating nc files we can use:
list_nc=`\ls nc/SELECTION_RGPS_S???_${YEAR}????h??_${YEAR}????h??.nc`
nbf=`echo ${list_nc} | wc -w`
echo " => ${nbf} files => ${nbf} batches!"


for ff in ${list_nc}; do

    fb=`basename ${ff}`
    echo
    echo " *** Doing file ${fb}"

    # Number of records inside netCDF file:
    Nr=`ncdump -h ${ff} | grep 'time\ =\ UNLIMITED' | cut -d'(' -f2 | cut -d' ' -f1`
    echo "   => it has ${Nr} time records"
    if [ ${Nr} -ne 2 ]; then
        echo "  ==> bad! Presently we expect 2 records!"
        exit
    fi

    
    for res in ${LIST_RES}; do

    
        CMD="${EXE} ${ff} ${res}"
        echo "    ==> will launch:"; echo "     ${CMD}"; echo
        ${CMD}
        echo
    done

done

#exit








#if [ "$1" = "d" ]; then
#    CMD="${EXE} ${FSI3IN} ${FNMM}"
#else
#CMD="${EXE}  ice_tracking.nc  ${FNMM} 0,72 330" ; # NEMO SEED at HSS=15



#fi

