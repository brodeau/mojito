#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=1

if not len(argv) in [4]:
    print('Usage: '+argv[0]+' <file_Q_mesh_N1.npz> <file_Q_mesh_N2.npz> <marker_size>')
    exit(0)
cf_Q1 = argv[1]
cf_Q2 = argv[2]
mrkrsz= int(argv[3])


# Reading the quad meshes in both npz files:
QUA1 = mjt.LoadClassPolygon( cf_Q1, ctype='Q' )
QUA2 = mjt.LoadClassPolygon( cf_Q2, ctype='Q' )


cdt1 = str(QUA1.date)
cdt2 = str(QUA2.date)
print('\n *** Dates for the two records:',cdt1,cdt2)

idt1 = clock2epoch(cdt1)
idt2 = clock2epoch(cdt2)

rdt = idt2 - idt1
print('      => `dt` = '+str(rdt)+' s, or '+str(rdt/(3600*24))+' days')


# Comprehensive name for npz and figs to save later on:
cf1, cf2 = path.basename(cf_Q1), path.basename(cf_Q2)
cnm_pref = split('_',cf1)[1]
ccdt1, ccdt2 = split('_',cdt1)[0], split('_',cdt2)[0]
cnm_pref = cnm_pref+'_'+ccdt1+'_'+ccdt2

print('\n *** Number of points in the two records:',QUA1.nP,QUA2.nP)
print('\n *** Number of quads in the two records:',QUA1.nQ,QUA2.nQ)




# The Quads we retain, i.e. those who exist in both snapshots:
vnm, vidx1, vidx2 = np.intersect1d( QUA1.QuadNames, QUA2.QuadNames, assume_unique=True, return_indices=True )
nQ = len(vnm) ; # also = len(vidx*) 
print('       => there are '+str(nQ)+' Quads common to the 2 records!\n')

# Coordinates of the 4 points of quadrangles for the 2 consecutive records:
zXY1 = QUA1.MeshPointXY[vidx1,:,:].copy() ; #* 1000.  ; # 1000 => from km to m
zXY2 = QUA2.MeshPointXY[vidx2,:,:].copy() ; #* 1000.

# Computation of partial derivative of velocity vector constructed from the 2 consecutive positions:
zX, zY, zU, zV, zdUdxy, zdVdxy = mjt.PDVfromPos( rdt, zXY1, zXY2, QUA1.area()[vidx1], QUA2.area()[vidx2],  iverbose=idebug )

# zX, zY => positions at center of time interval!
# zU, zV => velocities at center of time interval!

if idebug>0:
    # DEBUG: show velocities at each 4 vertices of each Quad:
    zzx = zX.flatten()
    zzy = zY.flatten()
    zzu = zU.flatten()
    zzv = zV.flatten()
    #
    mjt.ShowDeformation( zzx, zzy, zzu, cfig='./figs/'+cnm_pref+'_U4.png', cwhat='U4',
                         marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4 )
    mjt.ShowDeformation( zzx, zzy, zzv, cfig='./figs/'+cnm_pref+'_V4.png', cwhat='V4',
                         marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4 )
    #
    del zzx, zzy, zzu, zzv

# Coordinates of barycenter of Quads at center of time interval:
zXc = np.mean( zX[:,:], axis=1 )
zYc = np.mean( zY[:,:], axis=1 )

if idebug>1:
    # Velocities of barycenter of Quads at center of time interval:
    zUc = np.mean( zU[:,:], axis=1 )
    zVc = np.mean( zV[:,:], axis=1 )
    #
    mjt.ShowDeformation( zXc, zYc, zUc, cfig='./figs/'+cnm_pref+'_Uc.png', cwhat='Uc',
                         marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4 )
    mjt.ShowDeformation( zXc, zYc, zVc, cfig='./figs/'+cnm_pref+'_Vc.png', cwhat='Vc',
                         marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4 )
    #
    del zUc, zVc

# Divergence:
zdiv = mjt.DivPDV(zdUdxy, zdVdxy)

# Shear:
zshr = mjt.ShearPDV(zdUdxy, zdVdxy )

del zdUdxy, zdVdxy


# Saving data:
np.savez_compressed( './npz/DEFORMATIONS_'+cnm_pref+'.npz', Npoints=nQ, Xc=zXc, Yc=zYc, divergence=zdiv, shear=zshr )


# Some plots:
if not path.exists('./figs'): mkdir('./figs')

mjt.ShowDeformation( zXc, zYc, zdiv, cfig='./figs/'+cnm_pref+'_Divergence.png', cwhat='div',
                     marker_size=mrkrsz, pFmin=-1.e-5, pFmax=1.e-5, zoom=4 )
mjt.ShowDeformation( zXc, zYc, zshr, cfig='./figs/'+cnm_pref+'_Shear.png',      cwhat='shr',
                     marker_size=mrkrsz, pFmin=0.,      pFmax=0.8e-5,  zoom=4 )
