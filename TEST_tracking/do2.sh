#!/bin/sh

. ./conf.bash

EXE="${MOJITO_DIR}/generate_quad_mesh.py"


#if [ "$1" = "d" ]; then
#    CMD="${EXE} ${FSI3IN} ${FNMM}"
#else
CMD="${EXE}  ice_tracking.nc  ${FNMM} 0,72 330"
#fi

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo

${CMD}
