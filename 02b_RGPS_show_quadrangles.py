#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock
import lbrgps   as lbr

idebug=1

if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <file_Q_mesh.npz>')
    exit(0)
cf_npzQ = argv[1]

#cfroot = split('.npz',path.basename(cf_in))[0]

#cf_npzT = './npz/T-mesh_'+cfroot+'.npz'




print('\n *** Getting quad meshes in file '+cf_npzQ+' !')
    
dataQ       = np.load(cf_npzQ)  ; #, allow_pickle=True)
xCoor       = dataQ['pointCoordinates']
Quadrangles = dataQ['Quadrangles']

print('')

cf_fig = str.replace( path.basename(cf_npzQ), '.npz', '.png' )

# Show quadrangles on a map:
kk = lbr.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig=cf_fig,
                     QuadMesh=Quadrangles, lProj=False, zoom=5 )

