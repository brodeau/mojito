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
    print('Usage: '+argv[0]+' <file_Q_mesh_N1.npz> <file_Q_mesh_N2.npz>')
    exit(0)

cf_Q1 = argv[1]
cf_Q2 = argv[2]


print('')

#cf_fig = str.replace( path.basename(cf_Q1), '.npz', '.png' )

# Reading the quad meshes in the npz files:
QUA1 = lbr.LoadClassPolygon( cf_Q1, ctype='Q' )
QUA2 = lbr.LoadClassPolygon( cf_Q2, ctype='Q' )


print('\n nP:',QUA1.nP,QUA2.nP)
print('\n nQ:',QUA1.nQ,QUA2.nQ)

nP0 = QUA2.nP
nQ0 = QUA2.nQ

if nP0>QUA1.nP or nQ0>QUA1.nQ:
    print('ERROR: more points or quadrangles in second record/file!!! :()'); exit(0)


# We need indices (in QUA1 arrays) of points and quads that make it to second record/file
#for jQ in range(nQ0):
#    print(QUA1.PointNames[jQ], ' | ', QUA2.PointNames[jQ])

vnm1, vnm2 = QUA1.PointNames.copy(), QUA2.PointNames.copy()

# The Quads we retain, i.e. those who exist in the 2 snapshots:
vnm, vidx1, vidx2 = np.intersect1d( vnm1, vnm2, assume_unique=True, return_indices=True )
nQ = len(vnm) ; # also = len(vidx*) 

Q1coor, Q2coor = np.zeros((nQ,4,2)), np.zeros((nQ,4,2))
Q1coor[:,:,:] = QUA1.MeshPointXY[vidx1,:,:]
Q2coor[:,:,:] = QUA2.MeshPointXY[vidx2,:,:]


#print(xx)

print(len(vidx1),len(vidx2))

exit(0)

# Getting the angles:
zAngles  = QUA1.angles()
zLengths =  QUA1.lengths()

for jQ in range(nQ):
    print('\n *** Quad #'+str(jQ)+' => angles =',zAngles[jQ,:],' => lengths =',zLengths[jQ,:])




# Show the quads with only the points that define them:
#kk = lbr.ShowTQMesh( QUA1.PointXY[:,0], QUA1.PointXY[:,1], cfig=cf_fig, QuadMesh=QUA1.MeshPointIDs, lProj=False, zoom=izoom )
