#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
from glob import glob
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=1

cprefixIn='DEFORMATIONS_z' ; # Prefix of deformation files...

if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <directory_input_npz_files>')
    exit(0)
cd_in = argv[1]
#cf_Q2 = argv[2]
#mrkrsz= int(argv[3])

# Polpulating deformation files available:
listnpz = np.sort( glob(cd_in+'/'+cprefixIn+'*.npz') )
nbFiles = len(listnpz)
print('\n *** We found '+str(nbFiles)+' deformation files into '+cd_in+' !\n')


kStreamName = np.zeros(nbFiles, dtype='U4')
#kDate1      = np.zeros(nbFiles, dtype='U10')
kf = 0
for ff in listnpz:
    #print(ff)
    fb = path.basename(ff)
    vf = split('_',fb)
    #print(vf,'\n')
    kStreamName[kf] = vf[2]

    kf = kf+1

print('list of streams:', kStreamName[:])


### np.savez_compressed( './npz/DEFORMATIONS_'+cnm_pref+'.npz', Npoints=nQ, Xc=zXc, Yc=zYc, divergence=zdiv, shear=zshr )



exit(0)
