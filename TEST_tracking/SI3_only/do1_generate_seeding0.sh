#!/bin/bash

. ./conf.bash

EXE1="${MOJITO_DIR}/tracking/generate_idealized_seeding.py"
EXE2="python3 -u ${MOJITO_DIR}/tracking/si3_part_tracker.py"

echo $NDATE1
echo

fout="./nc/mojito_seeding_nemoTsi3_${NDATE1}_HSS${iHSS}.nc"

if [ ! -f ${fout} ]; then
    #${EXE1} ${LDATE1} ${FNMM},${iHSS}
    ${EXE1} ${LDATE1} ${FNMM},${iHSS} ${FSI3IN} 0
fi


# Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
CMD="${EXE2} ${FSI3IN} ${FNMM} ${fout}" ; # with nc file for init seed...
echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
