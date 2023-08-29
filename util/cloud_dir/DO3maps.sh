#!/bin/bash

. ./conf.bash

EXE="${MOJITO_DIR}/diags/plot_map_def.py"

dir_in="npz_maps"


VEXPS=( `echo ${EXPS} | sed -e s/','/' '/g` )
vexps=( ${VEXPS[*],,} ) ; # same in lower case
echo ${VEXPS[*]}

NO=`echo ${VEXPS[*]} | wc -w`
echo ${NO}
        
LEXPS=""
ji=0
while [ ${ji} -lt ${NO} ]; do

    echo ${ji}

    ff=`\ls  ${dir_in}/${vexps[${ji}]}/DEFORMATIONS_*${VEXPS[${ji}]}*_${BATCH}_dt120_SLCT10km.npz`
    echo ${ff}
    LEXPS+="${ff} "
    echo
    

    ji=$((ji+1))
done

echo " LEXPS = ${LEXPS} "; echo

for cf in "shear" "total" "divergence"; do
#for cf in "total"; do

    CMD="${EXE} ${cf} ${LEXPS}"
    
    echo
    echo " Will do:"
    echo ${CMD}
    ${CMD}
    echo

done
#exit
rsync -avP figs/maps/map_*.png ${DIR_EXPORT_CLOUD}/Maps/
