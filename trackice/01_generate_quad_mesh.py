#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#  INPUT DATA: a `npz` file created with `trackice/scripts/traj2npz.py` (conversion from CSV to NPZ)
#
#    L. Brodeau, August 2022
#
# TO DO: use `nemo_box = cp.nemo_hbox(CNEMO,CBOX)` !!!
#
#
#  ABOUT input `npz` file:
#   * Name: should be of the form `NANUK4_ICE-BBM00_6h_19960101_19961031(_xxx).npz`
#
##################################################################

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np
from re import split

from scipy.spatial import Delaunay

import lbrgps   as lbr

from climporn import chck4f

from netCDF4 import Dataset

from climporn import epoch2clock


idebug=1

l_plot = True ; # Create figures to see what we are doing...

l_work_with_dist = True ; # work with distance (x,y, Cartesian coordinates) rather than geographic coordinates (lon,lat)...

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

MinDistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]

# Selection of appropriate quadrangles:
rTang_min =  10. ; # minimum angle tolerable in a triangle [degree]
rTang_max = 120. ; # maximum angle tolerable in a triangle [degree]
#
rQang_min =  65.  ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 115.  ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.7 ; # value that `1 - abs(L/H)` should not overshoot!
rQarea_min =  0. ; # min area allowed for Quadrangle [km^2]

rQarea_max = 500. ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:1
#rQarea_max = 7000. ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:5

rzoom_fig = 10



if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    if len(argv) != 4:
        print('Usage: '+argv[0]+' <file_trj.npz> <LSM_file> <record to use (C)>')
        exit(0)

    cf_npz = argv[1]
    cf_lsm = argv[2]
    irec   = int(argv[3])

    chck4f(cf_npz)
    chck4f(cf_lsm)

    # Getting time info and time step from input npz file which is should look like NEMO output file:
    vv = split('-|_', path.basename(cf_npz))
    CCONF = vv[0]
    print('\n *** Original NEMO CONF = '+CCONF)

    # Need some calendar info:
    with np.load(cf_npz) as data:
        NrTraj = data['NrTraj']
        vtime  = data['time']
    #
    NbDays = int( (vtime[-1] - vtime[0]) / (3600.*24.) )
    cdt1 = epoch2clock(vtime[0] )
    cdt2 = epoch2clock(vtime[-1])
    #
    print('\n *** Trajectories contain '+str(NrTraj)+' records...')
    print('    *  start and End dates => '+cdt1+' -- '+cdt2)
    print('        ==> nb of days =', NbDays)
    #
    cdats = epoch2clock(vtime[irec])
    cdate = str.replace( epoch2clock(vtime[irec], precision='D'), '-', '')
    print('    *   will get data at record #'+str(irec)+'! => date =',cdats)

    vf = split('_|\.', path.basename(cf_npz))
    cfroot = vf[0]+'_'+vf[1]+'_'+vf[2]+'_'+vf[5]+'_'+vf[6]+'_'+cdate
    
    cf_npzT = './npz/T-mesh_'+cfroot+'.npz'
    cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'

    if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ)):

        # Have to build triangle and quadrangle mesh!

        print('\n *** We are going to build triangle and quad meshes!')

        #############################
        print('\n *** Reading into '+cf_npz+' !!!')
        with np.load(cf_npz) as data:
            xIDs    = data['IDs'][:,irec]
            xJIs    = data['JIs'][:,irec]
            xJJs    = data['JJs'][:,irec]
        
        (NbP,) = np.shape(xIDs)
        print('\n *** There are '+str(NbP)+' buoys at the begining...')

        # We need to load the NEMO's metric files to translate `jj,ji` to actual coordinates:
        print('\n *** Reading "'+CCONF+'" metrics in "'+cf_lsm+'" ...')
        with Dataset(cf_lsm) as id_lsm:
            nb_dim = len(id_lsm.variables['tmask'].dimensions)
            Ni = id_lsm.dimensions['x'].size
            Nj = id_lsm.dimensions['y'].size
            print('    --- the shape of the '+CCONF+' domain appears to be Ni, Nj =', Ni, Nj)
            #if not l_work_with_dist:
            xlon_t = id_lsm.variables['glamt'][0,:,:]
            xlat_t = id_lsm.variables['gphit'][0,:,:]
            xlon_u = id_lsm.variables['glamu'][0,:,:]
            xlat_v = id_lsm.variables['gphiv'][0,:,:]
            print('      done.')

        # Load `distance to coast` data:
        vlon_dist, vlat_dist, xdist = lbr.LoadDist2CoastNC( fdist2coast_nc )

        zlon, zlat = np.zeros(NbP), np.zeros(NbP)
        zmsk       = np.zeros(NbP, dtype=int)

        for jb in range(NbP):
            rjj, rji = xJJs[jb], xJIs[jb]

            jj, rj = int(rjj)-1, rjj%1.   ; # F2C !
            ji, ri = int(rji)-1, rji%1.   ; # F2C !
            # #fixme: I'm still not sure whether 0.5 means T point or U,V point !!! => find out !!!
            ####      => for now, assume T is at 0 and U is at 0.5 (that might be the opposite...)

            #  --  working with geographic coordinates rather than cartesian coordinates...
            if ri <= 0.5:
                # 0<=ri<=0.5 ==> then we must interpolate between T_i and U_i:
                rlon = 2.*(0.5-ri)*xlon_t[jj,ji] + 2.*ri*xlon_u[jj,ji]
            else:
                # 0.5<ri<1 ==> then we must interpolate between U_i and T_i+1:
                rlon = 2.*(1.-ri)*xlon_u[jj,ji] + 2*(ri-0.5)*xlon_t[jj,ji+1]

            if rj <= 0.5:
                # 0<=rj<=0.5 ==> then we must interpolate between T_j and V_j:
                rlat = 2.*(0.5-rj)*xlat_t[jj,ji] + 2.*rj*xlat_v[jj,ji]
            else:
                # 0.5<rj<1 ==> then we must interpolate between V_j and T_j+1:
                rlat = 2.*(1.-rj)*xlat_v[jj,ji] + 2*(rj-0.5)*xlat_t[jj+1,ji]

            rd_ini = lbr.Dist2Coast( rlon, rlat, vlon_dist, vlat_dist, xdist )

            if rd_ini > MinDistFromLand:
                zlon[jb] = rlon
                zlat[jb] = rlat
                zmsk[jb] = 1
                if idebug>1:
                    print(' jb =',jb,': rjj, rji =',rjj, rji)
                    print('      => zlat =', zlon[jb])
                    print('      => zlon =', zlat[jb])
                    print('')

        # Okay not we have to get rid of points that are masked (only viscinity to land at play so far...)
        zlon = np.ma.masked_where( zmsk==0, zlon)
        zlat = np.ma.masked_where( zmsk==0, zlat)
        zIDs = np.ma.masked_where( zmsk==0, xIDs)

        # Delete masked elements:
        zlon = np.ma.MaskedArray.compressed(zlon)
        zlat = np.ma.MaskedArray.compressed(zlat)
        zIDs = np.ma.MaskedArray.compressed(zIDs)

        # Update number of valid points:
        NbP = len(zlon)
        if len(zlat)!=NbP or len(zIDs)!=NbP:
            print('ERROR: len(zlat)!=NbP or len(zIDs)!=NbP !'); exit(0)
        print('\n *** We retain only '+str(NbP)+' buoys! (others were too close to land)')

        # Name for each point:
        vPnam = np.array( [ str(i) for i in zIDs ], dtype='U32' )
        print('vPnam =', vPnam)

        if l_work_with_dist:
            from cartopy.crs import PlateCarree, NorthPolarStereo
            crs_src = PlateCarree()
            crs_trg = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
            zx,zy,_ = crs_trg.transform_points(crs_src, zlon, zlat).T
            xCoor = np.array( [ zx, zy ] ).T / 1000. ; # to km
            del zx, zy
        else:
            xCoor = np.array( [ zlon, zlat ] ).T
        #
        del zlon, zlat, xJIs, xJJs, xIDs

        if idebug>0:
            for jc in range(NbP):
                print(' * Name: "'+vPnam[jc]+'": ID='+str(zIDs[jc])+', x_coor='+str(round(xCoor[jc,0],2))+', y_coor='+str(round(xCoor[jc,1],2)))
            print('')



        # Generating triangular meshes out of the cloud of points:
        TRI = Delaunay(xCoor)

        xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

        (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles

        xNeighborIDs = TRI.neighbors.copy() ;  # shape = (Nbt,3)

        print('\n *** We have '+str(NbT)+' triangles!')

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

        del xQpnts, xQcoor, xCoor

        if not path.exists('./npz'): mkdir('./npz')
        
        # Save the triangular mesh info:
        lbr.SaveClassPolygon( cf_npzT, TRIAS, ctype='T', date=cdats )

        # Save the quadrangular mesh info:
        lbr.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q', date=cdats )

        del TRIAS, QUADS

    #if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ))
    ############################################################


    if l_plot:
    
        # Reading the triangle and quad class objects in the npz files:
        TRI = lbr.LoadClassPolygon( cf_npzT, ctype='T' )
        QUA = lbr.LoadClassPolygon( cf_npzQ, ctype='Q' )


        #zA = QUA.area()
        #for jQ in range(len(zA)):
        #    print('Quad Area =', zA[jQ])
        #del zA

        
    
        if not path.exists('./figs'): mkdir('./figs')

        # Show triangles on a map:
        print('\n *** Launching Triangle plot!')
        kk = lbr.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/fig01_Mesh_Triangles_'+cfroot+'.png',
                             TriMesh=TRI.MeshPointIDs, lProj=(not l_work_with_dist), zoom=rzoom_fig)
    
        # Show triangles together with the quadrangles on a map:
        print('\n *** Launching Triangle+Quad plot!')
        kk = lbr.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/fig02_Mesh_Quadrangles_'+cfroot+'.png',
                             TriMesh=TRI.MeshPointIDs,
                             pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshPointIDs,
                             lProj=(not l_work_with_dist), zoom=rzoom_fig)
    
        ## Show only points composing the quadrangles:
        #kk = lbr.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfroot+'.png',
        #                     lProj=(not l_work_with_dist) )
    
        # Show only the quads with only the points that define them:
        print('\n *** Launching Quad-only plot!')
        kk = lbr.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfroot+'.png',
                             QuadMesh=QUA.MeshPointIDs, lProj=(not l_work_with_dist), zoom=rzoom_fig)

