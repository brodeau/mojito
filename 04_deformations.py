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

dt =  3600. * 3.   ;  # time interval between 2 records #fixme: use the real time !!!!



if not len(argv) in [3]:
    print('Usage: '+argv[0]+' <file_Q_mesh_N1.npz> <file_Q_mesh_N2.npz>')
    exit(0)

cf_Q1 = argv[1]
cf_Q2 = argv[2]

# Comprehensive name for npz and figs to save later on:
cf1, cf2 = path.basename(cf_Q1), path.basename(cf_Q2)
cnm_pref = split('_',cf1)[1]
cdt1, cdt2 = split('_',cf1)[4], split('_',cf2)[4]
cnm_pref = cnm_pref+'_'+cdt1+'-'+cdt2


# Reading the quad meshes in both npz files:
QUA1 = lbr.LoadClassPolygon( cf_Q1, ctype='Q' )
QUA2 = lbr.LoadClassPolygon( cf_Q2, ctype='Q' )

print('\n *** Number of points in the two records:',QUA1.nP,QUA2.nP)
print('\n *** Number of quads in the two records:',QUA1.nQ,QUA2.nQ)

# The Quads we retain, i.e. those who exist in both snapshots:
vnm, vidx1, vidx2 = np.intersect1d( QUA1.PointNames, QUA2.PointNames, assume_unique=True, return_indices=True )
nQ = len(vnm) ; # also = len(vidx*) 
print('       => there are '+str(nQ)+' Quads common to the 2 records!\n')

# Coordinates of the 4 points of quadrangles for the 2 consecutive records:
zXY1 = QUA1.MeshPointXY[vidx1,:,:].copy() ; #* 1000.  ; # 1000 => from km to m
zXY2 = QUA2.MeshPointXY[vidx2,:,:].copy() ; #* 1000.

# Computation of partial derivative of velocity vector constructed from the 2 consecutive positions:
zX, zY, zdUdxy, zdVdxy = lbr.PDVfromPos( dt, zXY1, zXY2, QUA1.area()[vidx1], QUA2.area()[vidx2],  iverbose=idebug )

# zX, zY => positions at center of time interval!

ztp1, ztp2 = np.zeros(nQ), np.zeros(nQ)

# Coordinates of barycenter of Quads at center of time interval:
zXc = np.mean( zX[:,:], axis=1 )
zYc = np.mean( zY[:,:], axis=1 )

# Divergence:
zdiv = np.zeros(nQ)
zdiv[:] = zdUdxy[:,0] + zdVdxy[:,1]

# Shear:
zshr = np.zeros(nQ)
ztp1[:] = zdUdxy[:,0] - zdVdxy[:,1]
ztp2[:] = zdUdxy[:,1] + zdVdxy[:,0]
zshr[:] = np.sqrt( ztp1*ztp1 + ztp2*ztp2 )

del ztp1, ztp2, zdUdxy, zdVdxy


# Saving data:
np.savez_compressed( './npz/DEFORMATIONS_'+cnm_pref+'.npz', Npoints=nQ, Xc=zXc, Yc=zYc, divergence=zdiv, shear=zshr )


# Some plots:
if not path.exists('./figs'): mkdir('./figs')

lbr.ShowDeformation( zXc, zYc, zdiv, cfig='./figs/'+cnm_pref+'_Divergence.png', cwhat='div', pFmin=-5.e-6, pFmax=5.e-6, zoom=4 )
lbr.ShowDeformation( zXc, zYc, zshr, cfig='./figs/'+cnm_pref+'_Shear.png',      cwhat='shr', pFmin=0.,     pFmax=1.e-5, zoom=4 )


