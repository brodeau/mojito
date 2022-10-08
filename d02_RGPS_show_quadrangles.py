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

if not len(argv) in [3]:
    print('Usage: '+argv[0]+' <file_Q_mesh.npz> <zoom>')
    exit(0)
cf_npzQ = argv[1]
izoom   = int(argv[2])

#cfroot = split('.npz',path.basename(cf_in))[0]

#cf_npzT = './npz/T-mesh_'+cfroot+'.npz'




print('\n *** Getting quad meshes in file '+cf_npzQ+' !')
    
dataQ = np.load(cf_npzQ)  ; #, allow_pickle=True)
XY    = dataQ['PointXY']
Mesh = dataQ['Mesh']

print('Shape pf `XY` =', np.shape(XY))
print('Shape pf `Mesh` =', np.shape(Mesh))


print('')

cf_fig = str.replace( path.basename(cf_npzQ), '.npz', '.png' )

# Show quadrangles on a map:
kk = lbr.ShowTQMesh( XY[:,0], XY[:,1], cfig=cf_fig, QuadMesh=Mesh, lProj=False, zoom=izoom )

