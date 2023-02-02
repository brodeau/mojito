#!/bin/bash

. ./conf.bash

#dt1=19970101
#dt2=19970104

EXE1="${MOJITO_DIR}/trackice/01_generate_quad_mesh.py"
EXE2="${MOJITO_DIR}/deformation.py"

fin=`basename ${FILIN} | cut -d '.' -f1`
cprf=`echo ${fin} | cut -d '_' -f1-3`_`echo ${fin} | cut -d '_' -f6-7`

#fQ1="./npz/Q-mesh_${cprf}_${dt1}.npz"
#fQ2="./npz/Q-mesh_${cprf}_${dt2}.npz"


#if [ ! -f ${fQ1} ] && [ ! -f ${fQ1} ]; then

Nrec=15
irec=0

while [ ${irec} -le ${Nrec} ]; do
    ir1=${irec}
    ir2=$((irec+1))
    exit;#lolo

    CMD="${EXE1} ${FILIN} ${MESH_MASK} 0,1 ${RESOL0}"

    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo
    ${CMD}

#fi
    

if [ ! -f ${fQ1} ] && [ ! -f ${fQ1} ]; then
    echo "ERROR: cannot find ${fQ1} and/or ${fQ2} !"
    exit
fi


echo; echo
echo " *** Computation of deformartions"
CMD="${EXE2} ${fQ1} ${fQ2} 0 ${MARKER_SIZE}" ; # time deviation has to be `0` because model output !!!
echo "  ==> ${CMD}"; echo
${CMD}
echo; echo
