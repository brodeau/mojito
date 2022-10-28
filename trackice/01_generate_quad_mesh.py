#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#  INPUT DATA: a `npz` file created with `trackice/scripts/traj2npz.py` (conversion from CSV to NPZ)
#
#    L. Brodeau, 2022
#
# #TODO: in ri,rj part make sure that ji or jj must not be updated (for lat after lon has been done) when one of these 2 are > 0.5
#                ex: ri = 0.9 ok keep ji,jj for lon stuff, but then should we not look at lat(ji+1,jj) for lat stuff ???
#
# #TODO: This script should be used for first record, but then, at time records later in time,
#        this script or another should follow the same exact quads found with 1st record, and not
#        start a Delaunay from scratch at this particular record !!!!
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

l_debug_plot = False

l_plot = True ; # Create figures to see what we are doing...

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
#rQarea_max = 200. ; # max area allowed for Quadrangle [km^2]
#
#rzoom_fig = 10


#rQarea_max = 500. ; rzoom_fig = 10 ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:1
#rQarea_max = 7000. ; rzoom_fig = 4  ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:5
#rQarea_max = 70000. ; rzoom_fig = 1  ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:5
rQarea_max = 18000. ; rzoom_fig = 2  ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:5







if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    if len(argv) != 4:
        print('Usage: '+argv[0]+' <file_trj.npz> <LSM_file> <records to use (C), comma-separated>')
        exit(0)

    cf_npz = argv[1]
    cf_lsm = argv[2]
    lirec = argv[3]

    vcrec = split(',',lirec)
    Nrec  = len(vcrec)
    vrec = [ int(vcrec[i]) for i in range(Nrec) ]
    del vcrec

    chck4f(cf_npz)
    chck4f(cf_lsm)


    #########################################################################################################

    # First record:
    irec = vrec[0]

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
        print('\n *** Reading for record # '+str(irec)+' into '+cf_npz+' !!!')
        with np.load(cf_npz) as data:
            xIDs    = data['IDs'][:,irec]
            xJIs    = data['JIs'][:,irec]
            xJJs    = data['JJs'][:,irec]

        (NbP,) = np.shape(xIDs)
        print('\n *** There are '+str(NbP)+' buoys at record #'+str(irec)+'...')

        # We need to load the NEMO's metric files to translate `jj,ji` to actual coordinates:
        print('\n *** Reading "'+CCONF+'" metrics in "'+cf_lsm+'" ...')
        with Dataset(cf_lsm) as id_lsm:
            nb_dim = len(id_lsm.variables['tmask'].dimensions)
            Ni = id_lsm.dimensions['x'].size
            Nj = id_lsm.dimensions['y'].size
            print('    --- the shape of the '+CCONF+' domain appears to be Ni, Nj =', Ni, Nj)
            xlon_t = id_lsm.variables['glamt'][0,:,:]
            xlat_t = id_lsm.variables['gphit'][0,:,:]
            xlon_u = id_lsm.variables['glamu'][0,:,:]
            xlat_v = id_lsm.variables['gphiv'][0,:,:]
            print('      done.')

        #if l_debug_plot:
        #    xCoor_dbg = np.zeros((NbP,2))
        #    vPids_dbg = np.zeros(NbP, dtype=int)

        NbP, xCoor, zPid, vPnam = lbr.rJIrJJtoCoord( xJJs, xJIs, xIDs, xlon_t, xlon_u, xlat_t, xlat_v,
                                                     rMinDistFromLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc )

        if idebug>0:
            for jc in range(NbP):
                print(' * #'+str(jc)+' => Name: "'+vPnam[jc]+'": ID='+str(zPid[jc])+', x_coor='+str(round(xCoor[jc,0],2))+', y_coor='+str(round(xCoor[jc,1],2)))
            print('')

        # Generating triangular meshes out of the cloud of points:
        TRI = Delaunay(xCoor)

        xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

        (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles

        xNeighborIDs = TRI.neighbors.copy() ;  # shape = (Nbt,3)

        print('\n *** We have '+str(NbT)+' triangles!')
        #
        # Conversion to the `Triangle` class:
        TRIAS = lbr.Triangle( xCoor, xTpnts, xNeighborIDs, zPid, vPnam )

        del xTpnts, xNeighborIDs, TRI

        # Merge triangles into quadrangles:
        xQcoor, vPids, xQpnts, vQnam = lbr.Tri2Quad( TRIAS, iverbose=idebug, anglRtri=(rTang_min,rTang_max),
                                                     ratioD=rdRatio_max, anglR=(rQang_min,rQang_max),
                                                     areaR=(rQarea_min,rQarea_max) )
        if len(xQpnts)<=0: exit(0)

        (NbQ,_) = np.shape(xQpnts)
        print('\n *** We have '+str(NbQ)+' quadrangles!')

        # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
        QUADS = lbr.Quadrangle( xQcoor, xQpnts, vPids, vQnam )

        del xQpnts, xQcoor, xCoor

        if not path.exists('./npz'): mkdir('./npz')

        # Save the triangular mesh info:
        lbr.SaveClassPolygon( cf_npzT, TRIAS, ctype='T', date=cdats )

        # Save the quadrangular mesh info:
        lbr.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q', date=cdats )

        #del TRIAS, QUADS
    else:

        # Reading the triangle and quad class objects in the npz files:
        TRIAS = lbr.LoadClassPolygon( cf_npzT, ctype='T' )
        QUADS = lbr.LoadClassPolygon( cf_npzQ, ctype='Q' )

        l_debug_plot = False ; # we don't have xCoor_dbg and vPids_dbg !

    #if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ))
    ############################################################


    if l_plot:
        if not path.exists('./figs'): mkdir('./figs')

        if l_debug_plot:
            # Show all initial points (out of TrackIce):
            print('\n *** Launching initial cloud point plot!')
            kk = lbr.ShowTQMesh( xCoor_dbg[:,0], xCoor_dbg[:,1], cfig='./figs/fig01a_cloud_points_'+cfroot+'.png',
                                 ppntIDs=vPids_dbg, lGeoCoor=False, zoom=rzoom_fig )

        # Show triangles on a map:
        print('\n *** Launching Triangle plot!')
        kk = lbr.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig='./figs/fig01_Mesh_Triangles_'+cfroot+'.png',
                             ppntIDs=TRIAS.PointIDs,
                             TriMesh=TRIAS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig)

        # Show triangles together with the quadrangles on a map:
        print('\n *** Launching Triangle+Quad plot!')
        kk = lbr.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig='./figs/fig02_Mesh_Quadrangles_'+cfroot+'.png',
                             ppntIDs=TRIAS.PointIDs, TriMesh=TRIAS.MeshVrtcPntIdx,
                             pX_Q=QUADS.PointXY[:,0], pY_Q=QUADS.PointXY[:,1], QuadMesh=QUADS.MeshVrtcPntIdx,
                             lGeoCoor=False, zoom=rzoom_fig)

        ## Show only points composing the quadrangles:
        #kk = lbr.ShowTQMesh( QUADS.PointXY[:,0], QUADS.PointXY[:,1], cfig='./figs/fig03a_Mesh_Points4Quadrangles_'+cfroot+'.png',
        #                     ppntIDs=QUADS.PointIDs, lGeoCoor=False, zoom=rzoom_fig )

        # Show only the quads with only the points that define them:
        print('\n *** Launching Quad-only plot!')
        kk = lbr.ShowTQMesh( QUADS.PointXY[:,0], QUADS.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfroot+'.png',
                             ppntIDs=QUADS.PointIDs,
                             QuadMesh=QUADS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig)

    ##############################################################################################################################















    print('\n\n#########################################################################')

    # 2nd record !!!
    ################
    irec = vrec[1]

    cdats = epoch2clock(vtime[irec])
    cdate = str.replace( epoch2clock(vtime[irec], precision='D'), '-', '')
    print('    *   will get data at record #'+str(irec)+'! => date =',cdats)

    vf = split('_|\.', path.basename(cf_npz))
    cfroot = vf[0]+'_'+vf[1]+'_'+vf[2]+'_'+vf[5]+'_'+vf[6]+'_'+cdate

    cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'

    if not path.exists(cf_npzQ):

        print('\n *** Reading for record # '+str(irec)+' into '+cf_npz+' !!!')
        with np.load(cf_npz) as data:
            xIDs    = data['IDs'][:,irec]
            # We only want to read what was used at former record:
            xIDs,vind,_ = np.intersect1d(xIDs, zPid, return_indices=True) ; # retain only values,indices of `xIDs` that exist in zPid
            xJIs    = data['JIs'][vind,irec]
            xJJs    = data['JJs'][vind,irec]

        # Here it's simpler we do not have to go through all the "point cloud" => "Delaunay" => "Triangles" => "Quads" process
        #  => we just rebuild the same quads out of the info we had from former records


        # Loop along the nQ quads found at former record
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nPprev = QUADS.nP
        nQprev = QUADS.nQ

        print(' *** At former record, we had '+str(nQprev)+' Quads relying on '+str(nPprev)+' points!')

        vPids = QUADS.PointIDs

        nQ = 0
        lst_PntIDs   = []
        lst_Qprev_ok = []

        for jQ in range(nQprev):

            l_mv_next = False

            print('\n --- Quad #'+str(jQ)+':')

            # Indices of the 4 points:
            v4p = QUADS.MeshVrtcPntIdx[jQ,:]
            vIDs = np.array( [ vPids[i] for i in v4p ], dtype=int )
            if idebug>0: print(vIDs[:])

            for jid in vIDs:
                if not jid in xIDs:
                    print(' Nuhh! Point with ID '+str(jid)+' does not exist anymore in this record')
                    print('  => this Quad is forgotten...')
                    l_mv_next = True
                    break

            if not l_mv_next:
                # Ok here we now that this Quad still exist in this new record:
                nQ = nQ + 1
                print('  => still exists at current record :)!')
                lst_Qprev_ok.append(jQ)
                lst_PntIDs.append(vIDs)

        vQindKeep = np.array( lst_Qprev_ok )
        del lst_Qprev_ok

        if idebug>0:
            pids_prev = QUADS.PointIDs

        # At this stage if no point has disapeared between the 2 records we have:
        #   => zids == QUADS.PointIDs

        zids = np.unique(np.array(lst_PntIDs))
        if idebug: print('  => IDs of points still involved:', zids)
        nP = len(zids)
        del vind, vIDs, lst_PntIDs

        # We have now `nQ` surviving Quads with indices: `vQindKeep`
        print('\n *** At present record, we have '+str(nQ)+' Quads / '+str(nQprev)+' that survived!')
        print('        and  '+str(nP)+' points  / '+str(nPprev)+' involved...')

        # We can now reduce xJJs, xJIs, xIDs accordingly:
        xIDs,vPindKeep,_ = np.intersect1d(xIDs, zids, return_indices=True) ; # retain only values,indices of `xIDs` that exist in `zids`
        xJIs        = xJIs[vPindKeep]
        xJJs        = xJJs[vPindKeep]

        
        NbP, xQcoor, vPids, vPnam = lbr.rJIrJJtoCoord( xJJs, xJIs, xIDs, xlon_t, xlon_u, xlat_t, xlat_v,  rMinDistFromLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc )



        #vi = QUADS.MeshVrtcPntIdx
        #vx = QUADS.MeshVrtcPntIDs()
        #print(vx[11,:],' name:',QUADS.QuadNames[11],vi[11,:])
        #print(vx[56,:],' name:',QUADS.QuadNames[56],vi[56,:])
        #exit(0)
        
        # We must update `vQindKeep` in case we suppressed points in `rJIrJJtoCoord()` (due to proximity to coast):
        if NbP < nP:
            vIDgone = np.setdiff1d( xIDs, vPids ) ; # keep the values of `xIDs` that are not in `vPids`
            print('\n *** IDs of points that were supressed due to coast proximity:\n',vIDgone[:])
            # We must suppress Quads that were relying on these points:
            lQrm = []
            for jid in vIDgone:
                ([jind],) = np.where(QUADS.PointIDs==jid) ; # we want point INDEX that gives this point ID !!!
                for jQ in vQindKeep:
                    if jind in QUADS.MeshVrtcPntIdx[jQ,:]: lQrm.append(jQ)
            if len(lQrm) != nP-NbP:
                print('ERROR: len(lQrm)!=NbP-len(xIDs) !',len(lQrm),NbP-len(xIDs)); exit(0)
            print('   ==> Quads to supress (ID, name):\n',[ (QUADS.QuaIDs[i], QUADS.QuadNames[i]) for i in lQrm ])

            vQindKeep = np.setdiff1d( vQindKeep, np.array(lQrm) ) ; # rm lQrm from vQindKeep

            # Now, the problem is that supression of a Quad can also lead to the vanishing of up to 3 other points!
            #zvPIDs = QUADS.MeshVrtcPntIDs()
            vPIDs_o = np.unique(zvPIDs[:        ,:])
            vPIDs_n = np.unique(zvPIDs[vQindKeep,:])
            #nP2s = len(vPIDs_o) - len(vPIDs_n)
            vPIDsRM = np.setdiff1d( vPIDs_o, vPIDs_n )
            vPIDsRM = np.setdiff1d( vPIDsRM, vIDgone)
            #if nP2s > nP-NbP:
            if len(vPIDsRM)>0:
                #print('   ==> we also need to supress '+str(nP2s-(nP-NbP))+' extra points that are not in use anymore!')
                print('   ==> we also need to supress '+str(len(vPIDsRM))+' extra points that are not in use anymore!')
            print('   ==> IDs of these points to supress =',vPIDsRM)


            _,vind,_ = np.intersect1d(vPIDs_o, vPIDs_n, return_indices=True) ; # retain only values,indices of `vPIDs_o` that exist in `vPIDs_n`
            
            print(vind,len(vind))
            
            #print( np.shape(vPidx_o),np.shape(vPidx_n) )

            print('shape(xQcoor) =',np.shape(xQcoor))
            xQcoor = xQcoor[vind,:]
            vPids  = vPids[vind]
            exit(0)
            
        
        QUADS_NEW = lbr.Quadrangle( xQcoor, QUADS.MeshVrtcPntIdx[vQindKeep], vPids, QUADS.QuadNames[vQindKeep] )

        # Save the quadrangular mesh info:
        lbr.SaveClassPolygon( cf_npzQ, QUADS_NEW, ctype='Q', date=cdats )

    else:

        QUADS_NEW = lbr.LoadClassPolygon( cf_npzQ, ctype='Q' )


    #if not path.exists(cf_npzQ)
    ############################

    if l_plot:

        print('\n *** Launching Quad-only plot!')
        kk = lbr.ShowTQMesh( QUADS_NEW.PointXY[:,0], QUADS_NEW.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfroot+'.png',
                             ppntIDs=QUADS_NEW.PointIDs,
                             QuadMesh=QUADS_NEW.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig)



