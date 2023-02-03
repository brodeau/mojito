#!/bin/bash

irec0=0
Nrecs=11

. ./conf.bash

EXE1="${MOJITO_DIR}/trackice/01_generate_quad_mesh.py"
EXE2="${MOJITO_DIR}/deformation.py"

fin=`basename ${FILIN} | cut -d '.' -f1`
cprf=`echo ${fin} | cut -d '_' -f1-3`_`echo ${fin} | cut -d '_' -f6-7`


echo ${cprf}

mkdir -p ./logs

ijob=0

irec=${irec0}
while [ ${irec} -lt ${Nrecs} ]; do
    
    irA=${irec}
    irB=$((irec+1))


    flogprf="${cprf}_recs_${irA}-${irB}"
    
    CMD="${EXE1} ${FILIN} ${MESH_MASK} ${irA},${irB} ${RESOL0}"
    
    echo ; echo " *** About to launch:"; echo "     ${CMD}"; echo
    
    ${CMD} 1>./logs/${flogprf}.out 2>./logs/${flogprf}.err &
    
    ijob=$((ijob+1))

    if [ $((ijob%NJPAR)) -eq 0 ]; then
        echo "Waiting! (ijob = ${ijob})...."
        wait
        echo; echo
    fi
    
    irec=$((irec+1))
    
done

exit

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
