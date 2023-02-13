#!/bin/bash

irec0=0

. ./conf.bash

EXE="${MOJITO_DIR}/trackice/01_generate_quad_mesh.py"

fin=`basename ${FILIN} | cut -d '.' -f1`
cprf=`echo ${fin} | cut -d '_' -f1-3`_`echo ${fin} | cut -d '_' -f6-7`


echo ${cprf}

mkdir -p ./logs

ijob=0

irec=${irec0}
while [ ${irec} -lt ${NbRec} ]; do
    
    irA=${irec}
    irB=$((irec+1))


    flogprf="${cprf}_recs_${irA}-${irB}"
    
    CMD="${EXE} ${FILIN} ${MESH_MASK} ${irA},${irB} ${RESOL0}"
    
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
echo



