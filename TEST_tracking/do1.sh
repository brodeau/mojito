#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/ice_part_tracker.py"


# 1/ populate the proper NC files to seed from:
dirin="${MOJITO_DIR}/TEST_brgps/nc"; echo ${dirin}

list_nc=`\ls ${dirin}/SELECTION_buoys_RGPS_*_${YEAR}????h??_${YEAR}????h??.nc`

nbf=`echo ${list_nc} | wc -w`

echo " *** We have ${nbf} files !"


for fnc in ${list_nc}; do

    fb=`basename ${fnc}`
    echo "   * File: ${fb} :"
    #cdr=`echo ${fb} | cut -d'_' -f5-8`
    #cdr=`echo ${cdr} | cut -d'.' -f1`
    #cdt1=`echo ${cdr} | cut -d'_' -f1`
    #cdt2=`echo ${cdr} | cut -d'_' -f2`
    #Ymmdd1=`echo ${cdt1} | cut -d'h' -f1`
    #hh1=`echo ${cdt1} | cut -d'h' -f2`
    #Ymmdd2=`echo ${cdt2} | cut -d'h' -f1`
    #hh2=`echo ${cdt2} | cut -d'h' -f2`    
    #echo " * Ymmdd1,  hh1 = ${Ymmdd1}, ${hh1}"
    #echo " * Ymmdd2,  hh2 = ${Ymmdd2}, ${hh2}"

    # Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
    CMD="${EXE} ${FSI3IN} ${FNMM} ${fnc}" ; # with nc file for init seed...
    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo
    ${CMD}



    echo; echo;
    echo "FIXME: ice_part_tracker.py has to be able to used the correct initial date!!! For now it is only 01/01 ???"
    exit
done

exit

if [ ! -f ${FNCSEED} ]; then
    echo " ERROR: file ${FNCSEED} is missing !"; exit
fi


if [ "$1" = "d" ]; then
    CMD="${EXE} ${FSI3IN} ${FNMM} ${DATE2}"
else
    CMD="${EXE} ${FSI3IN} ${FNMM} ${DATE2} ${FNCSEED}" ; # with nc file for init seed...
fi

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
