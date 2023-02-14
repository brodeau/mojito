#!/bin/sh

. ./conf.bash


if [ "$1" = "" ]; then
    echo "USAGE: $0 <file_pos_mojito.nc>"
    exit
fi

FIN="$1"
    



EXE="${MOJITO_DIR}/generate_quad_mesh.py"


#if [ "$1" = "d" ]; then
#    CMD="${EXE} ${FSI3IN} ${FNMM}"
#else
#CMD="${EXE}  ice_tracking.nc  ${FNMM} 0,72 330" ; # NEMO SEED at HSS=15

CMD="${EXE} ${FIN} 0,1 10" ; # RGPS seed

#fi

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo

${CMD}
