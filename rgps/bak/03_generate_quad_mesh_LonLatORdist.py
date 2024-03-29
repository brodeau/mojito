#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock
import mojito   as mjt

idebug=1

l_work_with_dist = True ; # work with distance (x,y, Cartesian coordinates) rather than geographic coordinates (lon,lat)...

# Selection of appropriate quadrangles:
rTang_min =  10. ; # minimum angle tolerable in a triangle [degree]
rTang_max = 120. ; # maximum angle tolerable in a triangle [degree]
#
rQang_min =  65.  ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 115.  ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.7 ; # value that `1 - abs(L/H)` should not overshoot!
rQarea_min =  70. ; # min area allowed for Quadrangle [km^2]
rQarea_max = 130. ; # max area allowed for Quadrangle [km^2]

rzoom_fig = 5



if __name__ == '__main__':

    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <SELECTION_buoys_RGPS_streamXXX_XXX.npz>')
        exit(0)
    cf_npz = argv[1]


    cc = '_gc'
    if l_work_with_dist: cc = '_cc'

    cfroot = str.replace( split('.npz',path.basename(cf_npz))[0] , 'SELECTION_buoys_RGPS_','' )

    cf_npzT = './npz/T-mesh_'+cfroot+'.npz'
    cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'

    if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ)):

        # Have to build triangle and quadrangle mesh!

        print('\n *** We are going to build triangle and quad meshes!')

        print('\n *** Reading into '+cf_npz+' !!!')
        with np.load(cf_npz) as data:
            it     = data['itime']
            cdate  = str( data['date'] )
            Nbuoys = data['Npoints']
            vids   = data['vids']
            if l_work_with_dist:
                vx = data['vx']
                vy = data['vy']
                if len(vids) != len(vx) or len(vids) != len(vy): print('ERROR Y11!') ; exit(0)
            else:
                vlon = data['vlon']
                vlat = data['vlat']
                if len(vids) != len(vlon) or len(vids) != len(vlat): print('ERROR Y12!') ; exit(0)

        NbP = len(vids) ; # number of points
        if NbP != Nbuoys: print('ERROR: NbP != Nbuoys !'); exit(0)

        if len(vx)!=NbP or len(vy)!=NbP:      print('ERROR Y13!') ; exit(0)
        if len(vids) != len(np.unique(vids)): print('ERROR Y14!') ; exit(0)

        print('\n *** Stream at '+epoch2clock(it)+' => '+str(NbP)+' points!')

        vIDs  = np.array( vids )
        del vids

        if l_work_with_dist:
            xCoor = np.array( [ [vx[i],vy[i]] for i in range(NbP) ] ) ; # original x,y cartesian cordinates of the RGPS data!
            #                                                               # => which is in Polar Stereographic projection, lon_0=-45, lat_ts=70
        else:
            xCoor = np.array( [ [vlon[i],vlat[i]] for i in range(NbP) ] ) ; # lon,lat projection used by Anton => applying reverse projection " "

        # Name for each point:
        vPnam = np.array( [ str(i) for i in vIDs ], dtype='U32' )

        if idebug>0:
            for jc in range(NbP):
                print(' * Name: "'+vPnam[jc]+'": ID='+str(vIDs[jc])+', x_coor='+str(round(xCoor[jc,0],2))+', y_coor='+str(round(xCoor[jc,1],2)))
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


        # Generating triangular meshes out of the cloud of points:
        TRI = Delaunay(xCoor)

        xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

        (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles

        xNeighborIDs = TRI.neighbors.copy() ;  # shape = (Nbt,3)

        print('\n *** We have '+str(NbT)+' triangles!')

        # Conversion to the `Triangle` class:
        TRIAS = mjt.Triangle( xCoor, xTpnts, xNeighborIDs, vIDs, vPnam ) ; #lolo

        del xTpnts, xNeighborIDs, TRI

        # Merge triangles into quadrangles:
        xQcoor, vPQids, xQpnts, vQnam = mjt.Tri2Quad( TRIAS, iverbose=idebug, anglRtri=(rTang_min,rTang_max),
                                                      ratioD=rdRatio_max, anglR=(rQang_min,rQang_max),
                                                      areaR=(rQarea_min,rQarea_max) )
        if len(xQpnts)<=0: exit(0)

        (NbQ,_) = np.shape(xQpnts)
        print('\n *** We have '+str(NbQ)+' quadrangles!')

        # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
        QUADS = mjt.Quadrangle( xQcoor, xQpnts, vPQids, vQnam, date=cdate )

        del xQpnts, xQcoor, xCoor

        # Save the triangular mesh info:
        mjt.SaveClassPolygon( cf_npzT, TRIAS, ctype='T' )

        # Save the quadrangular mesh info:
        mjt.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q' )

        del TRIAS, QUADS

    #if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ))
    ############################################################


    # Reading the triangle and quad class objects in the npz files:
    TRI = mjt.LoadClassPolygon( cf_npzT, ctype='T' )
    QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )


    if not path.exists('./figs'): mkdir('./figs')

    # Show triangles on a map:
    kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/fig01_Mesh_Triangles_'+cfroot+cc+'.png',
                         TriMesh=TRI.MeshVrtcPntIdx, lGeoCoor=(not l_work_with_dist), zoom=rzoom_fig)

    # Show triangles together with the quadrangles on a map:
    kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/fig02_Mesh_Quadrangles_'+cfroot+cc+'.png',
                         TriMesh=TRI.MeshVrtcPntIdx,
                         pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshVrtcPntIdx,
                         lGeoCoor=(not l_work_with_dist), zoom=rzoom_fig)

    ## Show only points composing the quadrangles:
    #kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfroot+cc+'.png',
    #                     lGeoCoor=(not l_work_with_dist) )

    # Show only the quads with only the points that define them:
    kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfroot+cc+'.png',
                         QuadMesh=QUA.MeshVrtcPntIdx, lGeoCoor=(not l_work_with_dist), zoom=rzoom_fig)

