#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/trackice/01_generate_quad_mesh.py"

ln -sf ${MESH_MASK} .

CMD="${EXE} ${FILIN} ${MESH_MASK} 0,1"

echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}

