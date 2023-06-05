#!/bin/bash

. ./conf.bash


if [ "$1" = "" ]; then
    echo "USAGE: $0 <Nb. batches>"
    exit
fi

Nb=${1}

jb=0

while [ ${jb} -le ${Nb} ]; do

    sb="S`printf "%03d" jb`"

    echo
    echo " *** Looking at batch # ${sb}"


done
