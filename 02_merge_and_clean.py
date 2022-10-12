#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################



from sys import argv, exit
from os import path
from glob import glob
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock
import lbrgps   as lbr

idebug=1



if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <dir_npz_selection>')
    exit(0)
cd_in = argv[1]

# Gathering all dates available from the files found:
vdates = []
listnpz = glob(cd_in+'/SELECTION_buoys_RGPS*.npz')
for ff in listnpz:        
    cdate = split( '_', split('.npz',path.basename(ff))[0] )[-1]
    print(' file, date =',ff,cdate)
    vdates.append(cdate)    

vdates = np.unique(vdates)
print('\n *** Dates available:',vdates[:],'\n')

# Now, for each date, going to merge the data from several files to one:
for cdate in vdates:
    print('    +++ date = ', cdate)
    listnpz = glob(cd_in+'/SELECTION_buoys_RGPS*_'+cdate+'.npz')
    nbf = len(listnpz)
    print('       '+str(nbf)+' files =>',listnpz)

    if nbf>1:
        print('         => will merge these '+str(nbf)+'!')
        kk = lbr.mergeNPZ( listnpz, cd_in+'/merges_selection_'+cdate+'.npz' )
        


exit(0)

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
    del vids
    
    if l_work_with_dist:
        xCoor = np.array( [ [vx[i]  ,vy[i]  ] for i in range(NbP) ] ) ; # original x,y cartesian cordinates of the RGPS data!
        #                                                               # => which is in Polar Stereographic projection, lon_0=-45, lat_ts=70
    else:
        xCoor = np.array( [ [vlon[i],vlat[i]] for i in range(NbP) ] ) ; # lon,lat projection used by Anton => applying reverse projection " "
    vnam  = np.array([ str(i) for i in vIDs ], dtype='U32')

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

    vPnam = np.array( [ str(i) for i in vIDs ], dtype='U32' )

    # Conversion to the `Triangle` class:
    TRIAS = lbr.Triangle( xCoor, xTpnts, xNeighborIDs, vPnam )

    del xTpnts, xNeighborIDs, TRI



    # Merge triangles into quadrangles:
    xQcoor, xQpnts, vQnam = lbr.Tri2Quad( TRIAS, iverbose=idebug, anglRtri=(rTang_min,rTang_max),
                                          ratioD=rdRatio_max, anglR=(rQang_min,rQang_max),
                                          areaR=(rQarea_min,rQarea_max) )
    if len(xQpnts)<=0: exit(0)

    (NbQ,_) = np.shape(xQpnts)
    print('\n *** We have '+str(NbQ)+' quadrangles!')

    # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
    QUADS = lbr.Quadrangle( xQcoor, xQpnts, vQnam )    

    del xQpnts, xQcoor

    del xCoor

    # Save the triangular mesh info:
    lbr.SaveClassPolygon( cf_npzT, TRIAS, ctype='T' )

    # Save the quadrangular mesh info:
    lbr.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q' )

#if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ))
############################################################


# Reading the triangle and quad class objects in the npz files:
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

