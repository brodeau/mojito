#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
from os import path
import numpy as nmp
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock
import lbrgps   as lbr

idebug=2

l_work_with_dist = True ; # work with distance (x,y, Cartesian coordinates) rather than geographic coordinates (lon,lat)...


if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <file_RGPS.npz>')
    exit(0)
cf_in = argv[1]

data = nmp.load(cf_in)

it   = data['itime']
vlon = data['vlon']
vlat = data['vlat']
vids = data['vids']

if len(vids) != len(vlon) or len(vids) != len(vlat):
    print('ERROR Y1!')
    exit(0)

crfig = split('.npz',path.basename(cf_in))[0]
    
ct = epoch2clock(it)

print('\n *** Stream at '+ct)

NbP = len(vlon) ; # number of points

print('\n *** We have '+str(NbP)+' points!')


vIDs  = nmp.array( vids )
xCoor = nmp.array([[vlon[i],vlat[i]] for i in range(NbP) ]             )
vnam  = nmp.array([ str(i) for i in vids ], dtype='U32')
#print(vnam)


if idebug>0:
    for jc in range(NbP):
        print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', lat='+str(round(xCoor[jc,1],2))+', lon='+str(round(xCoor[jc,0],2)))
    print('')


if l_work_with_dist:
    from cartopy import crs
    srs_src = crs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; #rgps
    srs_trg = crs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=60) ; # nextsim?

    zx,zy,_ = srs_trg.transform_points(srs_src, xCoor[:,0], xCoor[:,1]).T ; # km
    xCoor[:,0] = zx[:]
    xCoor[:,1] = zy[:]

    del zx, zy




    
# Generating triangular meshes out of the cloud of points:
TRI = Delaunay(xCoor)

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
kk = lbr.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='01_'+crfig+'.png',
                     pnames=vnam, TriMesh=xTriangles, lProj=False, izoom=7 )

# Merge triangles into quadrangles:
xQuads = lbr.Triangles2Quads( xTriangles, xNeighbors, xCoor, vnam,  iverbose=idebug )

if len(xQuads) <= 0: exit(0)

(NbQ,_) = nmp.shape(xQuads)
print('\n *** We have '+str(NbQ)+' quadrangles!')

# Show quadrangles on a map:
kk = lbr.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='02_'+crfig+'.png',
                     pnames=vnam, TriMesh=xTriangles, QuadMesh=xQuads, lProj=False, izoom=7 )


