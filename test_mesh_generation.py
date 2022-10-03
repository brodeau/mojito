#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
#from os import path
import numpy as nmp

from scipy.spatial import Delaunay

import lbrgps   as lbr

idebug=2

l_work_with_dist = True ; # work with distance (x,y, Cartesian coordinates) rather than geographic coordinates (lon,lat)...
l_cartopy = True
#l_cartopy = False

#l_work_with_dist = False ; # work with distance (x,y, Cartesian coordinates) rather than geographic coordinates (lon,lat)...

#y_gre, x_gre = 45.184369, 5.734251


xcities = [] ; ip = 0
#
xcities.append({"ID":ip,"city":"Paris"  ,  "lat":48.835334, "lon":2.353824 }); ip=ip+1
xcities.append({"ID":ip,"city":"Rome"   ,  "lat":41.89,     "lon":12.49    }); ip=ip+1
xcities.append({"ID":ip,"city":"Andorra",  "lat":42.506939, "lon":1.521247 }); ip=ip+1
xcities.append({"ID":ip,"city":"Athen"  ,  "lat":37.984149, "lon":23.727984}); ip=ip+1
xcities.append({"ID":ip,"city":"Belgrad",  "lat":44.817813, "lon":20.456897}); ip=ip+1
xcities.append({"ID":ip,"city":"Berlin" ,  "lat":52.517037, "lon":13.388860}); ip=ip+1
xcities.append({"ID":ip,"city":"Bern"   ,  "lat":46.948271, "lon":7.451451 }); ip=ip+1
xcities.append({"ID":ip,"city":"London" ,  "lat":51.510433, "lon":-0.129711}); ip=ip+1
xcities.append({"ID":ip,"city":"Esbjerg",  "lat":55.477434, "lon":8.468160 }); ip=ip+1
xcities.append({"ID":ip,"city":"Brest",    "lat":48.389657, "lon":-4.481700}); ip=ip+1
xcities.append({"ID":ip,"city":"Tunis",    "lat":36.802481, "lon":10.168440}); ip=ip+1
xcities.append({"ID":ip,"city":"Madrid",   "lat":40.414060, "lon":-3.699336}); ip=ip+1
xcities.append({"ID":ip,"city":"Alger",    "lat":36.732591, "lon": 3.101878}); ip=ip+1
xcities.append({"ID":ip,"city":"Dublin",   "lat":53.341825, "lon":-6.267302}); ip=ip+1
xcities.append({"ID":ip,"city":"Porto",    "lat":41.144559, "lon":-8.622812}); ip=ip+1
xcities.append({"ID":ip,"city":"Tanger",   "lat":35.771494, "lon":-5.836354}); ip=ip+1
xcities.append({"ID":ip,"city":"Bergen",   "lat":60.378183, "lon": 5.333694}); ip=ip+1
xcities.append({"ID":ip,"city":"Stockholm","lat":59.310631, "lon":18.067388}); ip=ip+1
xcities.append({"ID":ip,"city":"Ghardaia", "lat":32.483352, "lon": 3.681163}); ip=ip+1

NbP = len(xcities) ; # number of points
print('\n *** We have '+str(NbP)+' cities!')

# Conversion from dictionary to Numpy arrays:
vIDs  = nmp.array([ xcities[jc]["ID"]                      for jc in range(NbP) ], dtype=int  )
xCoor = nmp.array([[xcities[jc]["lon"],xcities[jc]["lat"]] for jc in range(NbP) ]             )
vnam  = nmp.array([ xcities[jc]["city"]                    for jc in range(NbP) ], dtype='U32')

if idebug>0:
    for jc in range(NbP):
        print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', lat='+str(round(xCoor[jc,1],2))+', lon='+str(round(xCoor[jc,0],2)))
    print('')



if 1==1:


    if l_work_with_dist:
        # Distance, aka cartesian coordinates, not degrees... => [km]
        x0, y0 = xCoor[:,0], xCoor[:,1]
        if l_cartopy:
            from cartopy.crs import PlateCarree, NorthPolarStereo, epsg
            crs_src = PlateCarree()
            #crs_trg = NorthPolarStereo(central_longitude=nmp.min(xCoor[:,0]), true_scale_latitude=45) ; # Europe!, x-axis starts at westernmost city!
            # Alternative:
            #crs_trg = epsg(4326) ; # LatLon with WGS84 datum used by GPS units and Google Earth
            crs_trg = epsg(3035)  ; # Europe!
            zx,zy,_ = crs_trg.transform_points(crs_src, x0, y0).T ; # to km!

        else:
            import pyproj as proj
            crs_src = proj.Proj(init='epsg:4326') # LatLon with WGS84 datum used by GPS units and Google Earth
            crs_trg = proj.Proj(init='epsg:3035') # Europe ?
            zx,zy   = proj.transform(crs_src, crs_trg, x0, y0)  # to km...


        xCoor[:,0],xCoor[:,1] = zx/1000., zy/1000.

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

    #print('class TRIAS:')
    #print(TRIAS.length,'\n')
    #print(TRIAS.ID,'\n')
    #print(TRIAS.pointIDs,'\n')
    #print(TRIAS.pointCoor,'\n')
    #print(TRIAS.neighbors,'\n')



    if idebug>1:
        for jT in range(NbT):
            vpl = TRIAS.pointIDs[jT,:] ;
            print(' Triangle #'+str(jT)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"')
            print('    => neighbor triangles are:',TRIAS.neighbors[jT,:],'\n')

    cc = '_gc'
    if l_work_with_dist: cc = '_cc'



    # Merge triangles into quadrangles:
    xQpnts, xQcoor = lbr.Triangles2Quads( TRIAS.pointIDs, TRIAS.neighbors, xCoor, vnam,  iverbose=idebug )
    if len(xQpnts) <= 0: exit(0)

    (NbQ,_) = nmp.shape(xQpnts)
    print('\n *** We have '+str(NbQ)+' quadrangles!')

    # Conversion to the `Quadrangle` class:
    QUADS = lbr.Quadrangle( NbQ, xQpnts, xQcoor )
    #print('class QUADS:')
    #print(QUADS.length,'\n')
    #print(QUADS.ID,'\n')
    #print(QUADS.pointIDs,'\n')
    #print(QUADS.pointCoor,'\n')
    #exit(0)
    del xQpnts, xQcoor

# Show triangles on a map:
kk = lbr.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='01_Mesh_Map_TRIangles_Europe'+cc+'.png',
                     pnames=vnam, TriMesh=TRIAS.pointIDs, lProj=(not l_work_with_dist))

# Show quadrangles on a map:
kk = lbr.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='02_Mesh_Map_Quadrangles_Europe'+cc+'.png',
                     pnames=vnam, TriMesh=TRIAS.pointIDs, QuadMesh=QUADS.pointIDs, lProj=(not l_work_with_dist) )


