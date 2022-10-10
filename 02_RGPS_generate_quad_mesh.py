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

l_work_with_dist = True ; # work with distance (x,y, Cartesian coordinates) rather than geographic coordinates (lon,lat)...
#l_cartopy = True

# Selection of appropriate quadrangles:
rTang_min =  10. ; # minimum angle tolerable in a triangle [degree]
rTang_max = 120. ; # maximum angle tolerable in a triangle [degree]
#
rQang_min =  65.  ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 115.  ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.5 ; # value that `1 - abs(L/H)` should not overshoot!
rQarea_min =  80. ; # min area allowed for Quadrangle [km^2]
rQarea_max = 120. ; # max area allowed for Quadrangle [km^2]


if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <file_RGPS.npz>')
    exit(0)
cf_in = argv[1]


cc = '_gc'
if l_work_with_dist: cc = '_cc'

cfroot = split('.npz',path.basename(cf_in))[0]

cf_npzT = './npz/T-mesh_'+cfroot+'.npz'
cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'

if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ)):

    # Have to build triangle and quadrangle mesh!

    print('\n *** We are going to build triangle and quad meshes!')

    data = np.load(cf_in)

    it   = data['itime']
    vids = data['vids']
    if l_work_with_dist:
        vx = data['vx']
        vy = data['vy']
        if len(vids) != len(vx) or len(vids) != len(vy): print('ERROR Y11!') ; exit(0)
    else:
        vlon = data['vlon']
        vlat = data['vlat']
        if len(vids) != len(vlon) or len(vids) != len(vlat): print('ERROR Y12!') ; exit(0)

    ct = epoch2clock(it)

    print('\n *** Stream at '+ct)

    NbP = len(vids) ; # number of points

    print('\n *** We have '+str(NbP)+' points!')


    vIDs  = np.array( vids )
    if l_work_with_dist:
        xCoor = np.array( [ [vx[i]  ,vy[i]  ] for i in range(NbP) ] ) ; # original x,y cartesian cordinates of the RGPS data!
        #                                                               # => which is in Polar Stereographic projection, lon_0=-45, lat_ts=70
    else:
        xCoor = np.array( [ [vlon[i],vlat[i]] for i in range(NbP) ] ) ; # lon,lat projection used by Anton => applying reverse projection " "
    vnam  = np.array([ str(i) for i in vids ], dtype='U32')

    if idebug>0:
        for jc in range(NbP):
            print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', x_coor='+str(round(xCoor[jc,0],2))+', y_coor='+str(round(xCoor[jc,1],2)))
        print('')


    #if l_work_with_dist:
    #    # Distance, aka cartesian coordinates, not degrees... => [km]
    #    x0, y0 = xCoor[:,0], xCoor[:,1]
    #    if l_cartopy:
    #        from cartopy.crs import PlateCarree, NorthPolarStereo ;#, epsg
    #        crs_src = PlateCarree()
    #        crs_trg = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
    #        zx,zy,_ = crs_trg.transform_points(crs_src, x0, y0).T
    #    else:
    #        print('FIX ME `pyproj`!'); exit(0)
    #        import pyproj as proj
    #        crs_src = proj.Proj(init='epsg:4326') # LatLon with WGS84 datum used by GPS units and Google Earth
    #        crs_trg = proj.Proj(init='epsg:3035') # Europe ?
    #        zx,zy   = proj.transform(crs_src, crs_trg, x0, y0)
    #
    #    xCoor[:,0],xCoor[:,1] = zx/1000., zy/1000. ; # to km...
    #    del x0, y0, zx, zy


    #AAAAAAA

    # Generating triangular meshes out of the cloud of points:
    TRI = Delaunay(xCoor)

    xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

    (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles

    xNeighborIDs = TRI.neighbors.copy() ;  # shape = (Nbt,3)

    print('\n *** We have '+str(NbT)+' triangles!')

    zTcoor = np.array([ [ xCoor[i,:] for i in xTpnts[jT,:] ] for jT in range(NbT) ])

    # Conversion to the `Triangle` class:
    TRIAS = lbr.Triangle( xTpnts, zTcoor, xNeighborIDs )

    del xTpnts, zTcoor, xNeighborIDs, TRI



    # Merge triangles into quadrangles:
    xQpnts, xQcoor = lbr.Tri2Quad( TRIAS, xCoor, vnam,  iverbose=idebug, anglRtri=(rTang_min,rTang_max),
                                   ratioD=rdRatio_max, anglR=(rQang_min,rQang_max), areaR=(rQarea_min,rQarea_max) )
    if len(xQpnts) <= 0: exit(0)

    (NbQ,_) = np.shape(xQpnts)
    print('\n *** We have '+str(NbQ)+' quadrangles!')

    # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
    QUADS = lbr.Quadrangle( lbr.TriPntIDs2QuaPntIDs(xQpnts), xQcoor )

    del xQpnts, xQcoor

    del xCoor

    # Save the triangular mesh info:
    lbr.SaveClassPolygon( cf_npzT, TRIAS, ctype='T' )
    
    # Save the quadrangular mesh info:
    lbr.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q' )

#if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ))
############################################################


# Reading the triangle and quad meshes in the npz files:
TRI = lbr.LoadClassPolygon( cf_npzT, ctype='T' )
QUA = lbr.LoadClassPolygon( cf_npzQ, ctype='Q' )


# Show triangles on a map:
kk = lbr.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='fig01_Mesh_Map_TRIangles_Europe'+cc+'.png',
                     TriMesh=TRI.MeshPointIDs, lProj=(not l_work_with_dist), zoom=5)

# Show triangles together with the quadrangles on a map:
kk = lbr.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='fig02_Mesh_Map_Quadrangles_Europe'+cc+'.png',
                     TriMesh=TRI.MeshPointIDs,
                     pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshPointIDs,
                     lProj=(not l_work_with_dist), zoom=5)

## Show only points composing the quadrangles:
#kk = lbr.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='fig03_Mesh_Map_Points4Quadrangles_Europe'+cc+'.png',
#                     lProj=(not l_work_with_dist) )

# Show only the quads with only the points that define them:
kk = lbr.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='fig03_Mesh_Map_Points4Quadrangles_Europe'+cc+'.png',
                     QuadMesh=QUA.MeshPointIDs, lProj=(not l_work_with_dist), zoom=5)

