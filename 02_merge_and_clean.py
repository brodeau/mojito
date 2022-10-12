#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################



from sys import argv, exit
from os import path
from glob import glob
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock
import lbrgps   as lbr

idebug=1



if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <dir_npz_selection>')
    exit(0)
cd_in = argv[1]

# Gathering all dates available from the files found:
vdates = []
listnpz = glob(cd_in+'/SELECTION_buoys_RGPS*.npz')
for ff in listnpz:        
    cdate = split( '_', split('.npz',path.basename(ff))[0] )[-1]
    print(' file, date =',ff,cdate)
    vdates.append(cdate)    

vdates = np.unique(vdates)
print('\n *** Dates available:',vdates[:],'\n')

# Now, for each date, going to merge the data from several files to one:
for cdate in vdates:
    print('    +++ date = ', cdate)
    listnpz = glob(cd_in+'/SELECTION_buoys_RGPS*_'+cdate+'.npz')
    nbf = len(listnpz)
    print('       '+str(nbf)+' files =>',listnpz)

    if nbf>1:
        print('         => will merge these '+str(nbf)+'!')
        kk = lbr.mergeNPZ( listnpz, cd_in+'/merged_selection_'+cdate+'.npz', iverbose=idebug )
