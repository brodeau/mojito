#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/diags/plot_map_def.py"

dir_in="npz_rgps_map"


VEXPS=( `echo ${EXPS} | sed -e s/','/' '/g` )

echo ${VEXPS[*]}

NO=`echo ${VEXPS[*]} | wc -w`
echo ${NO}
        

fbbm=""
fevp=""

frgps=`\ls ${dir_in}/rgps/DEFORMATIONS_RGPS_S003_dt120_1997*_10km.npz 2>/dev/null`
LEXPS="${frgps} "
if [ ${NO} -ge 2 ]; then
    fbbm=`\ls ${dir_in}/bbm2302/DEFORMATIONS_NEMO-SI3_NANUK4_${VEXPS[1]}_S003_dt120_1997*_10km.npz 2>/dev/null`
    LEXPS+="${fbbm} "
fi
if [ ${NO} -ge 3 ]; then
    fevp=`\ls ${dir_in}/evp2302/DEFORMATIONS_NEMO-SI3_NANUK4_${VEXPS[2]}_S003_dt120_1997*_10km.npz 2>/dev/null`
    LEXPS+="${fevp}"
fi


echo " LEXPS = ${LEXPS} "; echo

for cf in "shear" "divergence" "total"; do

    CMD="${EXE} ${cf} ${LEXPS}"
    
    echo
    echo " Will do:"
    echo ${CMD}
    ${CMD}
    echo

done

rsync -avP figs/maps/map_*.png ~/Nextcloud/Laurent/Work/2023/paper_bbmsi3/Maps/
