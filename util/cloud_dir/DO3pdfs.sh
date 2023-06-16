#!/bin/bash

. ./conf.bash

EXE="${HOME}/DEV/mojito/diags/plot_comp_PDFs.py"


VEXPS=( `echo ${EXPS} | sed -e s/','/' '/g`)
#VRESKM=( `echo ${SCALES} | sed -e s/','/' '/g`)
VRESKM=( 10 )
VDTBIN=( `echo ${BINS}   | sed -e s/','/' '/g`)

str_exps=`echo ${EXPS} | sed -e s/','/'_'/g`

echo ${VEXPS[*]}
echo ${VRESKM[*]}
echo ${VDTBIN[*]}

#exit
mkdir -p figs

ic=0
for RESKM in ${VRESKM[*]}; do

    DT_BINS_H=${VDTBIN[${ic}]}

    echo; echo " *** RESKM=${RESKM} ; DT_BINS_H=${DT_BINS_H} !"

    for cf in "shear" "total" "divergence" "convergence" "absDiv"; do

        list_ff=""
        
        for EXP in ${VEXPS[*]}; do

            exp=${EXP,,}

            echo "EXP, exp = ${EXP}, ${exp}"

            #if [ "${exp}" == "rgps" ]; then
            cdt=dt${DT_BINS_H}
            #else
            #    cdt=dt0
            #fi

            ff=`find -L ./${exp} -name "PDF_*${EXP}*_${cdt}_${RESKM}km_*_${cf}.npz" 2>/dev/null`


            echo " *** ff = ${ff}"

            nbw=`echo ${ff} | wc -w`
            if [ ${nbw} -ne 1 ]; then
                echo "PROBLEM: with ${EXP} ${RESKM} ${cf}"
                echo " => all files = ${ff}"
                exit
            fi
            
            list_ff+="${ff} "
            
        done
        echo " list_ff = ${list_ff} !"
        CMD="${EXE} ${list_ff}"

        echo
        echo " *** Wil launch: ${CMD}"

        echo
        ${CMD}
        echo
    done
    ic=`expr ${ic} + 1`
done

echo
echo " * Exporting figures to ${DIR_EXPORT_CLOUD}/Figures/PDFs !"

rsync -avP ./figs/PDFs/*.png ${DIR_EXPORT_CLOUD}/Figures/PDFs/



