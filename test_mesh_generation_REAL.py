#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
#from os import path
import numpy as nmp

from scipy.spatial import Delaunay

from climporn import epoch2clock
import lbrgps   as lbr

idebug=2

cf_in = 'npz/SELECTION_buoys_RGPS_stream001_1997-01-04.npz'

data = nmp.load(cf_in)

it   = data['itime']
vlon = data['vlon']
vlat = data['vlat']
vids = data['vids']

#print(it)
#print(vlon)
if len(vids) != len(vlon) or len(vids) != len(vlat):
    print('ERROR Y1!')
    exit(0)


ct = epoch2clock(it)

print('\n *** Stream at '+ct)

NbP = len(vlon) ; # number of points

print('\n *** We have '+str(NbP)+' points!')


vIDs  = nmp.array( vids )
Xcoor = nmp.array([[vlon[i],vlat[i]] for i in range(NbP) ]             )
vnam  = nmp.array([ str(i) for i in vids ], dtype='U32')
#print(vnam)


if idebug>0:
    for jc in range(NbP):
        print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', lat='+str(round(Xcoor[jc,1],2))+', lon='+str(round(Xcoor[jc,0],2)))
    print('')


# Generating triangular meshes out of the cloud of points:
TRI = Delaunay(Xcoor)

xTriangles = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

(NbT,_) = nmp.shape(xTriangles) ; # NbT => number of triangles

xNeighbors = TRI.neighbors.copy() ;  # shape = (Nbt,3)

print('\n *** We have '+str(NbT)+' triangles!')


if idebug>1:
    for jx in range(NbT):
        vpl = xTriangles[jx,:] ; # 3 point indices forming the triangle
        print(' Triangle #'+str(jx)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"')
        print('    => neighbor triangles are:',xNeighbors[jx,:],'\n')


# Show triangles on a map:
kk = lbr.ShowTMeshMap( Xcoor[:,0], Xcoor[:,1], xTriangles, cfig="01_Mesh_Map_TRIangles_Europe.png",
                       pnames=vnam, vextent=[10,350, 75,90] )
exit(0)

# Merge triangles into quadrangles:
xQuads = lbr.Triangles2Quads( xTriangles, xNeighbors, Xcoor, vnam,  iverbose=idebug )

(NbQ,_) = nmp.shape(xQuads)
print('\n *** We have '+str(NbQ)+' quadrangles!')

# Show quadrangles on a map:
kk = lbr.ShowQMeshMap( Xcoor[:,0], Xcoor[:,1], xQuads, cfig="02_Mesh_Map_Quadrangles_Europe.png", pnames=vnam, TriMesh=xTriangles )

