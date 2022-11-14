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
from netCDF4 import Dataset

from climporn import chck4f, epoch2clock
import mojito   as mjt


idebug=1

l_debug_plot = False

l_plot = True ; # Create figures to see what we are doing...
vrngX = [-2200.,1600.]
vrngY = [-1000.,1800.]

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
rQarea_max = 7000. ; rzoom_fig = 4  ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:5
#rQarea_max = 18000. ; rzoom_fig = 2  ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:10



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
    lstrec = argv[3]

    vcrec = split(',',lstrec)
    Nrec  = len(vcrec)
    vRec = np.array( [ int(vcrec[i]) for i in range(Nrec) ], dtype=int )
    del vcrec
    
    chck4f(cf_npz)
    chck4f(cf_lsm)

    if not path.exists('./npz'): mkdir('./npz')
    if l_plot:
        if not path.exists('./figs'): mkdir('./figs')

    #########################################################################################################

    # Getting time info and time step from input npz file which is should look like NEMO output file:
    vfi   = split('_|\.', path.basename(cf_npz))
    cfstr = vfi[0]+'_'+vfi[1]+'_'+vfi[2]+'_'+vfi[5]+'_'+vfi[6]
    CCONF = vfi[0]
    print('\n *** Original NEMO CONF = '+CCONF)
    
    # Need some calendar info:
    with np.load(cf_npz) as data:
        NbRec = data['NbRec']
        vtime = data['time']

    NbDays = int( (vtime[-1] - vtime[0]) / (3600.*24.) )
    cdt1 = epoch2clock(vtime[0] )
    cdt2 = epoch2clock(vtime[-1])

    print('\n *** Trajectories contain '+str(NbRec)+' records...')
    print('    *  start and End dates => '+cdt1+' -- '+cdt2)
    print('        ==> nb of days =', NbDays)

    
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


    # First explore how many buoys and how many records we have based on point IDs:
    with np.load(cf_npz) as data:
        xMSK = data['mask'][:,:]
        xIDs = data['IDs'][:,:]

    (NbP0,Nrtot) = np.shape(xIDs)
    
    # Disapearing buoys along the specified records?
    idxBok = np.arange(NbP0) ; # default: they are all fine!
    if np.any(xIDs[:,vRec]<1):
        (idxB0,idxR0) = np.where(xMSK==0)
        (idxB1,idxR1) = np.where(xIDs <1)
        if np.sum( np.abs(idxB0-idxB1) ) != 0 or np.sum( np.abs(idxR0-idxR1) ):
            print('ERROR: np.where(xIDs[:,vRec]<1) != np.where(xMSK==0) !'); exit(0)
        idx_vnshd = np.unique(idxB1)
        nbrm = len(idx_vnshd)
        print('\n *** The '+str(nbrm)+' buoys, with following IDs, do not make it along the '+str(Nrec)+' records specified:',xIDs[idx_vnshd,0])
        print('   ==> will be ignored!')
        (idxBok,) = np.where( xMSK[:,np.max(vRec)] == 1 ) ;  # those OK at the last record are good to keep !
        print('    ==>  updating number of buoys from '+str(NbP0)+' to '+str(len(idxBok))+'!')
        NbP0 = len(idxBok)
        del idxB0,idxR0,idxB1,idxR1
    
    if Nrec > Nrtot:
        print('ERROR: you want to work with more records than there is !!!',Nrec,Nrtot); exit(0)                    
    if np.max(vRec) > Nrtot-1:
        print('ERROR: max(vRec) > Nrtot-1 !', np.max(vRec), Nrtot-1 ); exit(0)

    del xIDs
    print('')
    
    # Array with the shape coresponding to the # of records we want to read:
    xIDs = np.zeros((NbP0,Nrec), dtype=int)
    xJIs = np.zeros((NbP0,Nrec), dtype=int)
    xJJs = np.zeros((NbP0,Nrec), dtype=int)
    mask = np.zeros( NbP0      , dtype=int) + 1  ; # Mask to for "deleted" points (to cancel)
    
    #############################
    with np.load(cf_npz) as data:
        jr = 0
        for jrec in vRec:
            print(' * Reading for record # '+str(jrec)+' into '+cf_npz+' !!!')
            xIDs[:,jr] = data['IDs'][idxBok,jrec]
            xJIs[:,jr] = data['JIs'][idxBok,jrec]
            xJJs[:,jr] = data['JJs'][idxBok,jrec]
            jr = jr+1
    
    zGC  = np.zeros((NbP0,2,Nrec))
    zXY  = np.zeros((NbP0,2,Nrec))            
    zPnm = np.array( [ str(i) for i in xIDs[:,0] ], dtype='U32' ) ; # Name for each point, based on 1st record...
    
    for jr in range(Nrec):
        print('\n   * Record #'+str(vRec[jr])+':')

        # Translate rji,rjj from TracIce to lon,lat and x,y:
        zGC[:,:,jr], zXY[:,:,jr] = mjt.rJIrJJtoCoord( xJJs[:,jr], xJIs[:,jr], xIDs[:,jr], xlon_t, xlon_u, xlat_t, xlat_v )
        
        # Get rid of points to close to land (shrinks arrays!):
        mask[:] = mjt.MaskCoastal( zGC[:,:,jr], mask=mask[:], rMinDistFromLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc )

    # How many points left after elimination of buoys that get too close to land (at any record):
    NbP = np.sum(mask)
    print('\n *** '+str(NbP)+' / '+str(NbP0)+' points survived the dist2coast test => ', str(NbP0-NbP)+' points to delete!')

    # We do not need to keep a record dimension for IDs, but first ensure there's no fuck-up...
    zPntIDs = xIDs[:,0].copy()
    for jr in range(Nrec):
        if np.any(xIDs[:,jr]-zPntIDs[:] != 0):
            print(' Fuck Up!!!! => `xIDs` not consistent along the records!!!' ); exit(0)
    del xIDs

    zPnm, zPntIDs, zGC, zXY = mjt.ShrinkArrays( mask, zPnm, zPntIDs, zGC, zXY )
    
    if idebug>1:
        for jr in range(Nrec):
            print('\n  DEBUG *** Record #'+str(vRec[jr])+':')
            for jc in range(NbP):                    
                print( ' * #'+str(jc)+' => Name: "'+zPnm[jc]+'": ID=',zPntIDs[jc],', X=',zXY[jc,0,jr],', Y=',zXY[jc,1,jr],
                       ', lon=',zGC[jc,0,jr],', lat=',zGC[jc,1,jr] )
                if str(zPntIDs[jc])!=zPnm[jc]:
                    print(' Fuck Up!!!! => zPntIDs[jc], zPnm[jc] =',zPntIDs[jc], zPnm[jc] ); exit(0)
        print('')


    for jr in range(Nrec):
        jrec = vRec[jr]
        
        print('\n\n *** QUAD-MESH GENERATION => record #'+str(jrec)+': record = '+str(jr))

        cdats  = epoch2clock(vtime[jrec])
        cdate  = str.replace( epoch2clock(vtime[jrec], precision='D'), '-', '')
        cfbase = cfstr+'_'+cdate                        
        print('    * which is original record '+str(jrec)+' => date =',cdats,'=>',cfbase)
        
        cf_npzQ = './npz/Q-mesh_'+cfbase+'.npz'

        
        if jr == 0:
            
            print('\n *** Delaunay triangulation for 1st record!')

            cf_npzT = './npz/T-mesh_'+cfbase+'.npz'
            
            # Generating triangular meshes out of the cloud of points for 1st record:
            lOK = False
            while not lOK:
            
                TRI = Delaunay(zXY[:,:,jr])
                xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*
                (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles
                xNeighborIDs = TRI.neighbors.copy() ;  # shape = (Nbt,3)

                # For some reasons, here, `xTpnts` can involve less points than the number of points
                # fed into Delaunay (zXY), if it is the case we have to shrink `zXY`, `zPntIDs`, `zPnm`
                # accordingly...
                PntIdxInUse = np.unique( xTpnts.flatten() )
                NbPnew      = len(PntIdxInUse)
                lOK = ( NbPnew == NbP )
                
                if not lOK:
                    if NbPnew > NbP:
                        print('ERROR: `NbPnew > NbP` !!!'); exit(0)
                    print('\n *** Need to cancel the '+str(NbP-NbPnew)+' points that vanished in Delaunay triangulation!')
                    mask  = np.zeros(  NbP, dtype=int )
                    zPntIdx0 = np.arange( NbP, dtype=int )
                    _,ind2keep,_ = np.intersect1d(zPntIdx0, PntIdxInUse, return_indices=True); # retain only indices of `zPntIdx0` that exist in `PntIdxInUse`
                    mask[ind2keep] = 1
                    if np.sum(mask) != NbPnew:
                        print('ERROR: `np.sum(mask) != NbPnew` !!!'); exit(0)
                    zPnm, zPntIDs, zGC, zXY = mjt.ShrinkArrays( mask, zPnm, zPntIDs, zGC, zXY )
                    del zPntIdx0, mask, ind2keep
                    NbP = NbPnew
                    print('      => done! Only '+str(NbP)+' points left in the game...')                    
    
                del PntIdxInUse, NbPnew

            ### while not lOK
            print('\n *** We have '+str(NbT)+' triangles!')
            
            # Conversion to the `Triangle` class:
            TRIAS = mjt.Triangle( zXY[:,:,jr], xTpnts, xNeighborIDs, zPntIDs, zPnm )
            del xTpnts, xNeighborIDs, TRI

            # Merge triangles into quadrangles:
            xQcoor, vPids, xQpnts, vQnam = mjt.Tri2Quad( TRIAS, iverbose=idebug, anglRtri=(rTang_min,rTang_max),
                                                         ratioD=rdRatio_max, anglR=(rQang_min,rQang_max),
                                                         areaR=(rQarea_min,rQarea_max) )
            if len(xQpnts)<=0: exit(0)

            (NbQ,_) = np.shape(xQpnts)
            print('\n *** We have '+str(NbQ)+' quadrangles!')

            # Save the triangular mesh info:
            mjt.SaveClassPolygon( cf_npzT, TRIAS, ctype='T', date=cdats )

            # To be used for other record, indices of Points to keep for Quads:
            _,ind2keep,_ = np.intersect1d(zPntIDs, vPids, return_indices=True); # retain only indices of `zPntIDs` that exist in `vPids`


            # Plots only for jrec==0:
            if l_plot:
                if l_debug_plot:
                    # Show all initial points (out of TrackIce):
                    print('\n *** Launching initial cloud point plot!')
                    kk = mjt.ShowTQMesh( zXY_dbg[:,0], zXY_dbg[:,1], cfig='./figs/fig01a_cloud_points_'+cfbase+'.png',
                                         ppntIDs=vPids_dbg, lGeoCoor=False, zoom=rzoom_fig,
                                         rangeX=vrngX, rangeY=vrngY )
                # Show triangles on a map:
                print('\n *** Launching Triangle plot!')
                kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig='./figs/fig01_Mesh_Triangles_'+cfbase+'.png',
                                     ppntIDs=TRIAS.PointIDs,
                                     TriMesh=TRIAS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig,
                                     rangeX=vrngX, rangeY=vrngY )


            #######################################################################################################
        else:
            #  => jr >= 1 !!!

            print('\n *** Quad tracking for jr=',jr)
            
            # Loop along the nQ quads found at former record
            nPprev = QUADS.nP
            nQprev = QUADS.nQ
    
            print(' *** At former record, we had '+str(nQprev)+' Quads relying on '+str(nPprev)+' points!')
    
            nQ = 0
            lst_PntIDs   = []
            for jQ in range(nQprev):
    
                l_Survived = True
    
                print('\n --- Quad #'+str(jQ)+':')
    
                # Indices of the 4 points:
                v4Pind = QUADS.MeshVrtcPntIdx[jQ,:]
                v4Pids = np.array( [ vPids[i] for i in v4Pind ], dtype=int )
    
                for jid in v4Pids:
                    if not jid in zPntIDs:
                        print(' Nuhh! Point with ID '+str(jid)+' does not exist anymore in this record')
                        print('  => this Quad is forgotten...')
                        l_Survived = False
                        break
    
                if l_Survived:
                    # Current Quad still exists in this record...
                    nQ = nQ + 1
                    print('     => still exists at current record :)!')
                    lst_PntIDs.append(v4Pids)
    
            zIDs = np.unique(np.array(lst_PntIDs))
            if idebug: print('  => IDs of points still involved:', zIDs)
            nP = len(zIDs)
            del v4Pind, v4Pids, lst_PntIDs, zIDs
    
            print('\n *** At present record, we have '+str(nQ)+' Quads / '+str(nQprev)+' that survived!')
            print('        and  '+str(nP)+' points  / '+str(nPprev)+' involved...')

            # * xQpnts remains the same! That's the whole point!!!
            # * xQcoor should be updated with the new point coordinates at this record:

            xQcoor[:,:] = zXY[ind2keep,:,jr]
            
        ### if jr == 0

        # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
        QUADS = mjt.Quadrangle( xQcoor, xQpnts, vPids, vQnam )

        # Save the quadrangular mesh info:
        mjt.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q', date=cdats )

        print('LOLO: `SaveClassPolygon` saved file',cf_npzQ,'with date =',cdats)
        

        if l_plot:
            # Show triangles together with the quadrangles on a map:
            print('\n *** Launching Triangle+Quad plot!')
            kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig='./figs/fig02_Mesh_Quadrangles_'+cfbase+'.png',
                                 ppntIDs=TRIAS.PointIDs, TriMesh=TRIAS.MeshVrtcPntIdx,
                                 pX_Q=QUADS.PointXY[:,0], pY_Q=QUADS.PointXY[:,1], QuadMesh=QUADS.MeshVrtcPntIdx,
                                 lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
            ## Show only points composing the quadrangles:
            #kk = mjt.ShowTQMesh( QUADS.PointXY[:,0], QUADS.PointXY[:,1], cfig='./figs/fig03a_Mesh_Points4Quadrangles_'+cfbase+'.png',
            #                     ppntIDs=QUADS.PointIDs, lGeoCoor=False, zoom=rzoom_fig )
            # Show only the quads with only the points that define them:
            print('\n *** Launching Quad-only plot!')
            kk = mjt.ShowTQMesh( QUADS.PointXY[:,0], QUADS.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                 ppntIDs=QUADS.PointIDs,
                                 QuadMesh=QUADS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig,
                                 rangeX=vrngX, rangeY=vrngY )
    
    ### for jr in range(Nrec):
    print('\n --- Over!\n')
