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

idebug=0
l_plot = True ; # Create figures to see what we are doing...

#l_box_restriction=True
l_box_restriction=False



if l_box_restriction:
    vrngX = [-400.,400.]
    vrngY = [-400.,400.]
else:
    vrngX = [-2200.,1600.]
    vrngY = [-1000.,1800.]

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

MinDistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]

# Selection of appropriate quadrangles:
#rTang_min =  10. ; # minimum angle tolerable in a triangle [degree]
#rTang_max = 120. ; # maximum angle tolerable in a triangle [degree]
rTang_min =   5. ; # minimum angle tolerable in a triangle [degree]
rTang_max = 160. ; # maximum angle tolerable in a triangle [degree]
#
#rQang_min =  65.  ; # minimum angle tolerable in a quadrangle [degree]
#rQang_max = 115.  ; # maximum angle tolerable in a quadrangle [degree]
rQang_min =  30.  ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 160.  ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.8 ; # value that `max(h1/h2,h2/h1)-1` should not overshoot! h1 being the "height" and "width" of the quadrangle

rzoom_fig = 4

if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    if len(argv) != 5:
        print('Usage: '+argv[0]+' <file_trj.npz> <LSM_file> <records to use (C), comma-separated> <basic_res_km>')
        exit(0)

    cf_npz = argv[1]
    cf_lsm = argv[2]
    lstrec = argv[3]
    resol0 = float(argv[4])

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


    if resol0>7.5 and resol0<12.5:
        rQarea_min =  60. ; # min area allowed for Quadrangle [km^2]
        rQarea_max = 200. ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:1
        #
    elif resol0>20 and resol0<30:
        rQarea_min = 350.
        rQarea_max = 950.
        #
    elif resol0>30 and resol0<40:
        rQarea_min =  800.
        rQarea_max = 1800.
        #
    elif resol0>50 and resol0<70:
        rQarea_min = 2000.
        rQarea_max = 5000.        
        #
    else:
        print('ERROR: we do not know what to do with resolution `resol0` =', resol0) ; exit(0)
        #
        #rQarea_max = 18000. ; rzoom_fig = 2  ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:10
        
        


    
    # Getting time info and time step from input npz file which is should look like NEMO output file:
    vfi   = split('_|\.', path.basename(cf_npz))
    cfstr = vfi[0]+'_'+vfi[1]+'_'+vfi[2]+'_'+vfi[5]+'_'+vfi[6]
    CCONF = vfi[0]
    print('\n *** Original NEMO CONF = '+CCONF)

    # Need some calendar info:
    with np.load(cf_npz) as data:
        NbRec = data['NbRec']
        vtime0 = data['time']

    NbDays = int( (vtime0[-1] - vtime0[0]) / (3600.*24.) )
    cdt1 = epoch2clock(vtime0[0] )
    cdt2 = epoch2clock(vtime0[-1])

    print('\n *** Trajectories contain '+str(NbRec)+' records...')
    print('    *  start and End dates => '+cdt1+' -- '+cdt2)
    print('        ==> nb of days =', NbDays)

    if np.any(vRec>=NbRec):
        print('ERROR: some of the specified records # are >= '+str(NbRec)+'  !'); exit(0)


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
    vtim = np.zeros(Nrec)
    mask = np.zeros( NbP0      , dtype=int) + 1  ; # Mask to for "deleted" points (to cancel)

    #############################
    with np.load(cf_npz) as data:
        jr = 0
        for jrec in vRec:
            print(' * Reading for record # '+str(jrec)+' into '+cf_npz+' !!!')
            xIDs[:,jr] = data['IDs'][idxBok,jrec]
            xJIs[:,jr] = data['JIs'][idxBok,jrec]
            xJJs[:,jr] = data['JJs'][idxBok,jrec]
            vtim[jr]   = data['time'][jrec]
            jr = jr+1

    zGC  = np.zeros((NbP0,2,Nrec))
    zXY  = np.zeros((NbP0,2,Nrec))
    ztim = np.zeros((NbP0,Nrec))
    zPnm = np.array( [ str(i) for i in xIDs[:,0] ], dtype='U32' ) ; # Name for each point, based on 1st record...

    for jp in range(NbP0):
        ztim[jp,:] = vtim[:]
    del vtim
            
    for jr in range(Nrec):
        print('\n   * Record #'+str(vRec[jr])+':')

        # Translate rji,rjj from TracIce to lon,lat and x,y:
        zGC[:,:,jr], zXY[:,:,jr] = mjt.rJIrJJtoCoord( xJJs[:,jr], xJIs[:,jr], xIDs[:,jr], xlon_t, xlon_u, xlat_t, xlat_v )

        # Get rid of points to close to land (shrinks arrays!):
        mask[:] = mjt.MaskCoastal( zGC[:,:,jr], mask=mask[:], rMinDistLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc )

    # How many points left after elimination of buoys that get too close to land (at any record):
    NbP1 = np.sum(mask)
    print('\n *** '+str(NbP1)+' / '+str(NbP0)+' points survived the dist2coast test => ', str(NbP0-NbP1)+' points to delete!')

    NbP = NbP1
    if l_box_restriction:
        # Restriction to a smaller box for debugging purposes:
        idxban = np.where( zGC[:,1,0] < 85. ) ; # keep only points north of 85.degN
        mask[idxban] = 0
        NbP = np.sum(mask)
        print('\n *** '+str(NbP)+' / '+str(NbP1)+' points survived the box test => ', str(NbP1-NbP)+' points to delete!')

    
    # We do not need to keep a record dimension for IDs, but first ensure there's no fuck-up...
    zPntIDs = xIDs[:,0].copy()
    for jr in range(Nrec):
        if np.any(xIDs[:,jr]-zPntIDs[:] != 0):
            print(' Fuck Up!!!! => `xIDs` not consistent along the records!!!' ); exit(0)
    del xIDs

    zPnm, zPntIDs, zGC, zXY, ztim = mjt.ShrinkArrays( mask, zPnm, zPntIDs, zGC, zXY, ztim )

    if idebug>1:
        for jr in range(Nrec):
            print('\n  DEBUG *** Record #'+str(vRec[jr])+':')
            for jc in range(NbP):
                print( ' * #'+str(jc)+' => Name: "'+zPnm[jc]+'": ID=',zPntIDs[jc],', X=',zXY[jc,0,jr],', Y=',zXY[jc,1,jr],
                       ', lon=',zGC[jc,0,jr],', lat=',zGC[jc,1,jr], ', time=',epoch2clock(ztim[jc,jr]) )
                if str(zPntIDs[jc])!=zPnm[jc]:
                    print(' Fuck Up!!!! => zPntIDs[jc], zPnm[jc] =',zPntIDs[jc], zPnm[jc] ); exit(0)
        print('')


    
    cdate0  = str.replace( epoch2clock(vtime0[vRec[0]], precision='D'), '-', '')
    
    for jr in range(Nrec):
        jrec = vRec[jr]

        print('\n\n *** QUAD-MESH GENERATION => record #'+str(jrec)+': record = '+str(jr))

        cdats  = epoch2clock(vtime0[jrec])
        cdate  = str.replace( epoch2clock(vtime0[jrec], precision='D'), '-', '')

        cfbase = cfstr+'_'+cdate0+'t0_'+cdate
        print('    * which is original record '+str(jrec)+' => date =',cdats,'=>',cfbase)

        cf_npzQ = './npz/Q-mesh_'+cfbase+'.npz'
        #print(' cf_npzQ =', cf_npzQ); exit(0)
        
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
                    zPnm, zPntIDs, zGC, zXY, ztim = mjt.ShrinkArrays( mask, zPnm, zPntIDs, zGC, zXY, ztim )
                    del zPntIdx0, mask, ind2keep
                    NbP = NbPnew
                    print('      => done! Only '+str(NbP)+' points left in the game...')

                del PntIdxInUse, NbPnew

            ### while not lOK
            print('\n *** We have '+str(NbT)+' triangles!')

            # Conversion to the `Triangle` class:
            TRIAS = mjt.Triangle( zXY[:,:,jr], xTpnts, xNeighborIDs, zPntIDs, ztim[:,jr], zPnm )
            del xTpnts, xNeighborIDs, TRI

            # Merge triangles into quadrangles:
            xQcoor, vPids, vTime, xQpnts, vQnam = mjt.Tri2Quad( TRIAS, anglRtri=(rTang_min,rTang_max),
                                                                ratioD=rdRatio_max, anglR=(rQang_min,rQang_max),
                                                                areaR=(rQarea_min,rQarea_max), idbglev=idebug )
            if len(xQpnts)<=0: exit(0)

            (NbQ,_) = np.shape(xQpnts)
            print('\n *** We have '+str(NbQ)+' quadrangles!')

            # Save the triangular mesh info:
            mjt.SaveClassPolygon( cf_npzT, TRIAS, ctype='T' )

            # To be used for other record, indices of Points to keep for Quads:
            _,ind2keep,_ = np.intersect1d(zPntIDs, vPids, return_indices=True); # retain only indices of `zPntIDs` that exist in `vPids`


            # Plots only for jrec==0:
            if l_plot:
                if idebug>1:
                    # Show all initial points (out of TrackIce):
                    print('\n *** Launching initial cloud point plot!')
                    kk = mjt.ShowTQMesh( zXY_dbg[:,0], zXY_dbg[:,1], cfig='./figs/fig01a_cloud_points_'+cfbase+'.png',
                                         ppntIDs=vPids_dbg, lGeoCoor=False, zoom=rzoom_fig,
                                         rangeX=vrngX, rangeY=vrngY )
                #
                if idebug>0:
                    # Show triangles on a map:
                    print('\n *** Launching Triangle plot!')
                    kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig='./figs/fig01_Mesh_Triangles_'+cfbase+'.png',
                                         ppntIDs=TRIAS.PointIDs,
                                         TriMesh=TRIAS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig,
                                         rangeX=vrngX, rangeY=vrngY )
                    
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS0 = mjt.Quadrangle( xQcoor, xQpnts, vPids, vTime, vQnam, date=cdats )

            # Save the quadrangular mesh info:
            mjt.SaveClassPolygon( cf_npzQ, QUADS0, ctype='Q' )


            #######################################################################################################
        else:
            
            #  => jr >= 1 !!!

            print('\n *** Quad recycling for jr=',jr)

            # Recycling Quads found at 1st record (QUADS0): lilo
            xQcoor, vTime, xQpnts, vPids, vQnam, vQIDs  = mjt.RecycleQuads( zXY[:,:,jr], ztim[:,jr], zPntIDs, QUADS0,  iverbose=idebug )
            
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS = mjt.Quadrangle( xQcoor, xQpnts, vPids, vTime, vQnam, vQIDs=vQIDs, date=cdats )

            # Save the quadrangular mesh info:
            mjt.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q' )

            del QUADS

        ### if jr == 0


        if l_plot:

            TRI = mjt.LoadClassPolygon( cf_npzT, ctype='T' )
            QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )


            if idebug>0:
                # Show triangles together with the quadrangles on a map:
                print('\n *** Launching Triangle+Quad plot!')
                kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/fig02_Mesh_Quadrangles_'+cfbase+'.png',
                                     ppntIDs=TRI.PointIDs, TriMesh=TRI.MeshVrtcPntIdx,
                                     pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshVrtcPntIdx,
                                     lGeoCoor=False, zoom=rzoom_fig,
                                     rangeX=vrngX, rangeY=vrngY )
                
            ## Show only points composing the quadrangles:
            #kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/fig03a_Mesh_Points4Quadrangles_'+cfbase+'.png',
            #                     ppntIDs=QUA.PointIDs, lGeoCoor=False, zoom=rzoom_fig )
            # Show only the quads with only the points that define them:
            print('\n *** Launching Quad-only plot!')
            kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/fig03_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                 ppntIDs=QUA.PointIDs,
                                 QuadMesh=QUA.MeshVrtcPntIdx, qIDs=QUA.QuadIDs, lGeoCoor=False, zoom=rzoom_fig,
                                 rangeX=vrngX, rangeY=vrngY )

    ### for jr in range(Nrec):
    print('\n --- Over!\n')
