#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
from os import path, mkdir
import numpy as np

from scipy.spatial import Delaunay

import lbrgps   as lbr

idebug=3

l_work_with_dist = True ; # work with distance (x,y, Cartesian coordinates) rather than geographic coordinates (lon,lat)...
#l_cartopy = True
l_cartopy = False

#y_gre, x_gre = 45.184369, 5.734251

# Selection of appropriate quadrangles:
rTang_min =  10. ; # minimum angle tolerable in a triangle [degree]  # #fixme: should be linked to rdRatio_max somehow???
rTang_max = 120. ; # maximum angle tolerable in a triangle [degree]
#
rQang_min =  65.  ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 115.  ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.5 ; # value that `1 - abs(L/H)` should not overshoot!
rQarea_min = 0.   ; # min area allowed for Quadrangle [km^2]
rQarea_max = 8.e5 ; # max area allowed for Quadrangle [km^2]


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
vIDs  = np.array([ xcities[jc]["ID"]                      for jc in range(NbP) ], dtype=int  )
xCoor = np.array([[xcities[jc]["lon"],xcities[jc]["lat"]] for jc in range(NbP) ]             )
vnam  = np.array([ xcities[jc]["city"]                    for jc in range(NbP) ], dtype='U32')

if idebug>0:
    for jc in range(NbP):
        print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', lat='+str(round(xCoor[jc,1],2))+', lon='+str(round(xCoor[jc,0],2)))
    print('')

cc = '_gc'
if l_work_with_dist: cc = '_cc'

cf_npzT = './npz/T-mesh_Europe.npz'
cf_npzQ = './npz/Q-mesh_Europe.npz'

if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ)):

    if not path.exists('./npz'): mkdir('./npz')

    # Have to build triangle and quadrangle mesh!

    cAu = 'degrees^2'
    if l_work_with_dist:
        cAu = 'km^2'
        # Distance, aka cartesian coordinates, not degrees... => [km]
        x0, y0 = xCoor[:,0], xCoor[:,1]
        if l_cartopy:
            from cartopy.crs import PlateCarree, NorthPolarStereo, epsg
            crs_src = PlateCarree()
            #crs_trg = NorthPolarStereo(central_longitude=np.min(xCoor[:,0]), true_scale_latitude=45) ; # Europe!, x-axis starts at westernmost city!
            # Alternative:
            #crs_trg = epsg(4326) ; # LatLon with WGS84 datum used by GPS units and Google Earth
            crs_trg = epsg(3035)  ; # Europe!
            zx,zy,_ = crs_trg.transform_points(crs_src, x0, y0).T ; # to km!

        else:
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

    (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles

    xNeighbors = TRI.neighbors.copy() ;  # shape = (Nbt,3)

    print('\n *** We have '+str(NbT)+' triangles!')

    zTcoor = np.array([ [ xCoor[i,:] for i in xTpnts[jT,:] ] for jT in range(NbT) ])

    # Conversion to the `Triangle` class:
    TRIAS = lbr.Triangle( NbT, xTpnts, zTcoor, xNeighbors )

    del xTpnts, zTcoor, xNeighbors, TRI

    if idebug>2:
        # SOME DEBUG TESTS TO SEE IF EVERYTHING WORKS FINE IN THE TRIANGLE CLASS:
        zlengths = TRIAS.lengths()
        zangles = TRIAS.angles()
        zarea   = TRIAS.area()

        for jT in range(TRIAS.nT):
            vpl = TRIAS.MeshPointIDs[jT,:] ; # IDs of the 3 points composing triangle
            print(' Triangle #'+str(jT)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"')
            print('    => neighbor triangles are:',TRIAS.neighbors[jT,:])
            print('    =>  lengths =',zlengths[jT,:])
            print('    =>  angles =',zangles[jT,:])
            print('    =>  area   =',round(zarea[jT],1),cAu)
            print('')
        del zlengths, zangles, zarea

    lbr.SavePolygon( cf_npzT, TRIAS, ctype='T' )


    # Merge triangles into quadrangles:
    xQpnts, xQcoor = lbr.Tri2Quad( TRIAS, xCoor, vnam,  iverbose=idebug, anglRtri=(rTang_min,rTang_max),
                                   ratioD=rdRatio_max, anglR=(rQang_min,rQang_max), areaR=(rQarea_min,rQarea_max) )
    if len(xQpnts) <= 0: exit(0)

    (NbQ,_) = np.shape(xQpnts)
    print('\n *** We have '+str(NbQ)+' quadrangles!')

    # Conversion to the `Quadrangle` class:
    QUADS = lbr.Quadrangle( NbQ, xQpnts, xQcoor )

    del xQpnts, xQcoor

    if idebug>2:
        # SOME DEBUG TESTS TO SEE IF EVERYTHING WORKS FINE IN THE QUADRANGLE CLASS:
        zlengths = QUADS.lengths()
        zangles  = QUADS.angles()
        zarea    = QUADS.area()
        zXYcloud = QUADS.PointXY
        print('\n  *** DEBUG summary for Quadrangle class:')
        print('  ***************************************')
        print('    => number of Quadrangles: '+str(QUADS.nQ))
        print('    ==> number of points involved: '+str(QUADS.nP))
        #print('    => Cloud of points (size ='+str(len(zXYcloud[:,0]))+') =')
        #print( zXYcloud[:,0] )
        #print( zXYcloud[:,1],'\n' )
        print('')
        for jQ in range(QUADS.nQ):
            vpl     = QUADS.MeshPointIDs[jQ,:] ; # IDs of the 4 points composing the quadrangle
            i1=vpl[0]; i2=vpl[1]; i3=vpl[2]; i4=vpl[3]
            vpr     = QUADS.MeshPointTriIDs[jQ,:] ; # IDs of the 4 points composing the quadrangle BUT in the reduced cloud of points
            print(' Quadrangle #'+str(jQ)+': ', vpl[:],'aka "'+vnam[i1]+' - '+vnam[i2]+' - '+vnam[i3]+' - '+vnam[i4]+'"')
            print('    =>  Pt. IDs =',vpl,' (with IDs from initial triangle points)')
            print('    =>  Pt. IDs =',vpr,' (with IDs from remaining points of quads)')
            print('    =>  coord.  =', round(xCoor[i1,0],0),round(xCoor[i1,1],0),round(xCoor[i2,0],0),round(xCoor[i2,1],0),
                                       round(xCoor[i3,0],0),round(xCoor[i3,1],0),round(xCoor[i4,0],0),round(xCoor[i4,1],0)  )
            print('    =>  lengths =',zlengths[jQ,:])
            print('    =>  angles =',zangles[jQ,:])
            print('    =>  area   =',round(zarea[jQ],1),cAu)
            print('')

        del zlengths, zangles, zarea
        #exit(0);#lolo

    del xCoor

    # Save the quadrangular mesh info:
    lbr.SavePolygon( cf_npzQ, QUADS, ctype='Q' )    
    #np.savez( cf_npzQ, PointXY=QUADS.PointXY, Mesh=QUADS.MeshPointIDs, MeshTriPntIDs=QUADS.MeshPointTriIDs,
    #          Lengths=QUADS.lengths(), Angles=QUADS.angles(), Areas=QUADS.area(), names=vnam )
    #print('\n *** "'+cf_npzQ+'" written!')

    # For plot to come:
    #Triangles   = TRIAS.MeshPointIDs
    #xyT         = TRIAS.PointXY

    #Quadrangles = QUADS.MeshPointIDs
    #xyQ         = QUADS.PointXY
    #QuadsRstrct = QUADS.MeshPointTriIDs
    
########################################################



print('\n\n *** Reading the triangle and quad meshes in the npz files...')

dataT = np.load(cf_npzT, allow_pickle=True)
xyT      = dataT['PointXY']
#vnam       = dataT['names']
Triangles  = dataT['Mesh']

dataQ = np.load(cf_npzQ, allow_pickle=True)
xyQ         = dataQ['PointXY']
Quadrangles = dataQ['Mesh']
QuadsRstrct = dataQ['MeshTriPntIDs']
print('')



# Show triangles on a map:
kk = lbr.ShowTQMesh( xyT[:,0], xyT[:,1], cfig='fig01_Mesh_Map_TRIangles_Europe'+cc+'.png',
                     pnames=vnam, TriMesh=Triangles, lProj=(not l_work_with_dist))

# Show triangles together with the quadrangles on a map:
kk = lbr.ShowTQMesh( xyT[:,0], xyT[:,1], cfig='fig02_Mesh_Map_Quadrangles_Europe'+cc+'.png',
                     pnames=vnam, TriMesh=Triangles, QuadMesh=Quadrangles, lProj=(not l_work_with_dist) )

## Show only points composing the quadrangles:
#kk = lbr.ShowTQMesh( xyQ[:,0], xyQ[:,1], cfig='fig03_Mesh_Map_Points4Quadrangles_Europe'+cc+'.png',
#                     lProj=(not l_work_with_dist) )

# Show only the quads with only the points that define them:
kk = lbr.ShowTQMesh( xyQ[:,0], xyQ[:,1], cfig='fig03_Mesh_Map_Points4Quadrangles_Europe'+cc+'.png',
                     QuadMesh=QuadsRstrct, lProj=(not l_work_with_dist) )

