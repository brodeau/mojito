#!/bin/bash

. ./conf.bash

EXE1="${MOJITO_DIR}/tracking/generate_idealized_seeding.py"
EXE2="python3 -u ${MOJITO_DIR}/tracking/si3_part_tracker.py"

YYYY=`echo ${DATE1} | cut -c1-4`
MM=`echo ${DATE1} | cut -c5-6`
DD=`echo ${DATE1} | cut -c7-8`


NDATE1="${YYYY}-${MM}-${DD}"
echo $NDATE1

LDATE1="${NDATE1}_00:00:00"
echo

fout="./nc/mojito_seeding_nemo_Tpoint_${NDATE1}_HSS${iHSS}.nc"

if [ ! -f ${fout} ]; then
    ${EXE1} ${LDATE1} ${FNMM},${iHSS}
fi




# Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
CMD="${EXE2} ${FSI3IN} ${FNMM} ${fout}" ; # with nc file for init seed...
echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
