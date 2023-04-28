#!/bin/bash

. ./conf.bash

str_exps=`echo ${EXPS} | sed -e s/','/'_'/g`

fout="scaling_TotDefRate_${str_exps}.pdf"

EXE="${HOME}/DEV/mojito/mkscaling.py"

#CMD1="inkscape --export-text-to-path --export-type=pdf"


CMD="${EXE} `pwd` ${EXPS} ${SCALES}"

echo; echo " *** Wild launch: ${CMD}"; echo
${CMD}
echo


