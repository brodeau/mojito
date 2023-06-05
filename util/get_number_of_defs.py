#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
from glob import glob
import numpy as np
from re import split
#
import mojito   as mjt
from mojito import config as cfg


if __name__ == '__main__':

    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <DEFORMATION_xxx.npz>')
        exit(0)    

    cf_in  = argv[1]

    with np.load(cf_in) as data:
        #rdate = int( data['time'] )
        nPnts =      data['Npoints']
        #if kf==0:
        #    corigin = str(data['origin'])
        #    reskm   = int(data['reskm_nmnl'])
        #fb = path.basename(ff)
        #vf = split('_',fb)
    #
    print(nPnts)
    exit(0)

