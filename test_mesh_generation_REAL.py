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
l_cartopy = True


if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <file_RGPS.npz>')
    exit(0)
cf_in = argv[1]


cfroot = split('.npz',path.basename(cf_in))[0]

cf_npzT = './npz/T-mesh_'+cfroot+'.npz'
cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'


if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ)):

    # Have to build triangle and quadrangle mesh!

    print('\n *** We are going to build triangle and quad meshes!')

    data = nmp.load(cf_in)

    it   = data['itime']
    vlon = data['vlon']
    vlat = data['vlat']
    vids = data['vids']

    if len(vids) != len(vlon) or len(vids) != len(vlat):
        print('ERROR Y1!')
        exit(0)

    ct = epoch2clock(it)

    print('\n *** Stream at '+ct)

    NbP = len(vlon) ; # number of points

    print('\n *** We have '+str(NbP)+' points!')


    vIDs  = nmp.array( vids )
    xCoor = nmp.array([[vlon[i],vlat[i]] for i in range(NbP) ]             )
    vnam  = nmp.array([ str(i) for i in vids ], dtype='U32')

    if idebug>0:
        for jc in range(NbP):
            print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', lat='+str(round(xCoor[jc,1],2))+', lon='+str(round(xCoor[jc,0],2)))
        print('')


    if l_work_with_dist:
        # Distance, aka cartesian coordinates, not degrees... => [km]
        x0, y0 = xCoor[:,0], xCoor[:,1]
        if l_cartopy:
            from cartopy.crs import PlateCarree, NorthPolarStereo ;#, epsg
            crs_src = PlateCarree()
            crs_trg = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
            zx,zy,_ = crs_trg.transform_points(crs_src, x0, y0).T
        else:
            print('FIX ME `pyproj`!'); exit(0)
            import pyproj as proj
            crs_src = proj.Proj(init='epsg:4326') # LatLon with WGS84 datum used by GPS units and Google Earth
            crs_trg = proj.Proj(init='epsg:3035') # Europe ?
            zx,zy   = proj.transform(crs_src, crs_trg, x0, y0)

    xCoor[:,0],xCoor[:,1] = zx/1000., zy/1000. ; # to km...
    del x0, y0, zx, zy

    #AAAAAAA

    # Generating triangular meshes out of the cloud of points:
    TRI = Delaunay(xCoor)

    xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

    (NbT,_) = nmp.shape(xTpnts) ; # NbT => number of triangles

    xNeighbors = TRI.neighbors.copy() ;  # shape = (Nbt,3)

    print('\n *** We have '+str(NbT)+' triangles!')

    zTcoor = nmp.array([ [ xCoor[i,:] for i in xTpnts[jT,:] ] for jT in range(NbT) ])

    # Conversion to the `Triangle` class:
    TRIAS = lbr.Triangle( NbT, xTpnts, zTcoor, xNeighbors )

    del xTpnts, zTcoor, xNeighbors, TRI

    if idebug>1:
        for jT in range(NbT):
            vpl = TRIAS.TriPointIDs[jT,:] ;
            print(' Triangle #'+str(jT)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"')
            print('    => neighbor triangles are:',TRIAS.neighbors[jT,:],'\n')

    cc = '_gc'
    if l_work_with_dist: cc = '_cc'

    # Merge triangles into quadrangles:
    xQpnts, xQcoor = lbr.Triangles2Quads( TRIAS.TriPointIDs, TRIAS.neighbors, xCoor, vnam,  iverbose=idebug )
    if len(xQpnts) <= 0: exit(0)

    (NbQ,_) = nmp.shape(xQpnts)
    print('\n *** We have '+str(NbQ)+' quadrangles!')

    # Conversion to the `Quadrangle` class:
    QUADS = lbr.Quadrangle( NbQ, xQpnts, xQcoor )
    #print('class QUADS:')
    #print(QUADS.length,'\n')
    #print(QUADS.ID,'\n')
    #print(QUADS.QuaPointIDs,'\n')
    #print(QUADS.pointCoor,'\n')
    #exit(0)
    del xQpnts, xQcoor

    # Save the triangular mesh info:
    nmp.savez( cf_npzT, pointCoordinates=xCoor, Triangles=TRIAS.TriPointIDs, names=vnam )
    print('\n *** "'+cf_npzT+'" written!')
    
    # Save the quadrangular mesh info:
    nmp.savez( cf_npzQ, pointCoordinates=xCoor, Quadrangles=QUADS.QuaPointIDs, names=vnam )
    print('\n *** "'+cf_npzQ+'" written!')

    # For plot to come:
    Triangles   = TRIAS.TriPointIDs
    Quadrangles = QUADS.QuaPointIDs

    
else:

    print('\n *** We are going to READ triangle and quad meshes in the npz files...')
    
    dataT = nmp.load(cf_npzT, allow_pickle=True)
    xCoor      = dataT['pointCoordinates']
    #vnam       = dataT['names']
    Triangles  = dataT['Triangles']
    
    dataQ = nmp.load(cf_npzQ, allow_pickle=True)
    Quadrangles = dataQ['Quadrangles']

    print('')


# Show triangles on a map:
kk = lbr.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='01_'+cfroot+'.png',
                     TriMesh=Triangles, lProj=(not l_work_with_dist), zoom=3)

# Show quadrangles on a map:
kk = lbr.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='02_'+cfroot+'.png',
                     TriMesh=Triangles, QuadMesh=Quadrangles, lProj=(not l_work_with_dist), zoom=3 )

