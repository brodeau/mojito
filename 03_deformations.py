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


print('')

#cf_fig = str.replace( path.basename(cf_Q1), '.npz', '.png' )

# Reading the quad meshes in the npz files:
QUA1 = lbr.LoadClassPolygon( cf_Q1, ctype='Q' )
QUA2 = lbr.LoadClassPolygon( cf_Q2, ctype='Q' )


print('\n *** nP:',QUA1.nP,QUA2.nP)
#print('\n nQ:',QUA1.nQ,QUA2.nQ)

nP0 = QUA2.nP
nQ0 = QUA2.nQ

if nP0>QUA1.nP or nQ0>QUA1.nQ:
    print('ERROR: more points or quadrangles in second record/file!!! :()'); exit(0)


# The Quads we retain, i.e. those who exist in the 2 snapshots:
vnm, vidx1, vidx2 = np.intersect1d( QUA1.PointNames, QUA2.PointNames, assume_unique=True, return_indices=True )
#print(len(vidx1),len(vidx2))
nQ = len(vnm) ; # also = len(vidx*) 
print('\n *** There are '+str(nQ)+' Quads common to the 2 records!')

# Coordinates of the 4 points of quadrangles:
#zXY1 , zXY2 = np.zeros((nQ,4,2)) - 999. , np.zeros((nQ,4,2)) - 999.
zXY1 = QUA1.MeshPointXY[vidx1,:,:].copy() ; #* 1000.  ; # 1000 => from km to m
zXY2 = QUA2.MeshPointXY[vidx2,:,:].copy() ; #* 1000.
# Same, but at center of time interval:
zX , zY = np.zeros((nQ,4)) , np.zeros((nQ,4))
zX = 0.5*( zXY1[:,:,0] + zXY2[:,:,0] )
zY = 0.5*( zXY1[:,:,1] + zXY2[:,:,1] )

# Area of quadrangles:
#zA1 , zA2 = np.zeros(nQ) - 999. , np.zeros(nQ) - 999.
#zA1[:] = QUA1.area()[vidx1]
#zA2[:] = QUA2.area()[vidx2]
# Area of quadrangles at center of time interval:
zA = np.zeros(nQ) - 999.
zA[:] = 0.5*( QUA1.area()[vidx1] + QUA2.area()[vidx2] ) ; #* 1.e6 ; # 1.e6 => from km^2 to m^2

# Velocities at center of time interval:
zU = np.array( [ zXY2[:,k,0] - zXY1[:,k,0] for k in range(4) ] ).T / dt ; # 1000 because X,Y in km !!!
zV = np.array( [ zXY2[:,k,1] - zXY1[:,k,1] for k in range(4) ] ).T / dt ; # 1000 because X,Y in km !!!
#zuu , zvv = np.zeros((nQ,4)) - 999. , np.zeros((nQ,4)) - 999.
#for jQ in range(nQ):
#    zuu[jQ,:] = 1000. * np.array( [ zXY2[jQ,k,0] - zXY1[jQ,k,0] for k in range(4) ] ) / dt
#    zvv[jQ,:] = 1000. * np.array( [ zXY2[jQ,k,1] - zXY1[jQ,k,1] for k in range(4) ] ) / dt
del zXY1, zXY2


if idebug>0:
    for jQ in range(0,nQ,100):
        print('  areas =',np.round(zA[jQ],3),'km^2, U =',np.round(zU[jQ,:],5),'m/s, V =',np.round(zV[jQ,:],5),'m/s')


# Partial derivatives:
#  --- the fact that units for coordinates was km and for area km^2 has no importance because it cancels,
#      we are looking to something in [s-1]
zdUdx , zdUdy = np.zeros(nQ) , np.zeros(nQ)
zdVdx , zdVdy = np.zeros(nQ) , np.zeros(nQ)
for jQ in range(nQ):
    zdUdx[jQ] =  np.sum( np.array([ (zU[jQ,(k+1)%4] + zU[jQ,k])*(zY[jQ,(k+1)%4] - zY[jQ,k]) for k in range(4) ]) ) / (2*zA[jQ])
    zdUdy[jQ] = -np.sum( np.array([ (zU[jQ,(k+1)%4] + zU[jQ,k])*(zX[jQ,(k+1)%4] - zX[jQ,k]) for k in range(4) ]) ) / (2*zA[jQ])
    zdVdx[jQ] =  np.sum( np.array([ (zV[jQ,(k+1)%4] + zV[jQ,k])*(zY[jQ,(k+1)%4] - zY[jQ,k]) for k in range(4) ]) ) / (2*zA[jQ])
    zdVdy[jQ] = -np.sum( np.array([ (zV[jQ,(k+1)%4] + zV[jQ,k])*(zX[jQ,(k+1)%4] - zX[jQ,k]) for k in range(4) ]) ) / (2*zA[jQ])

if idebug>0:
    for jQ in range(0,nQ,100):
        print('  dU/dx =',zdUdx[jQ],'1/s, dU/dy =',zdUdy[jQ],'1/s')
        print('  dV/dx =',zdVdx[jQ],'1/s, dV/dy =',zdVdy[jQ],'1/s\n')
    




        
# Show the quads with only the points that define them:
#kk = lbr.ShowTQMesh( QUA1.PointXY[:,0], QUA1.PointXY[:,1], cfig=cf_fig, QuadMesh=QUA1.MeshPointIDs, lProj=False, zoom=izoom )
