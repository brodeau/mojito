#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#  INPUT DATA: a `nc` file created with `ncSaveCloudBuoys()` of `mojito/ncio.py` !
#
#    L. Brodeau, 2023
#
##################################################################

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np
from re import split
from scipy.spatial import Delaunay

import mojito   as mjt

idebug = 0
iplot  = 1 ; # Create figures to see what we are doing...

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
#
rdev_scale = 0.25 ; # how much can we deviate from the specified scale to accept or reject a quadrangle
#                   # =>  (reskm*(1-rdev_scale))**2  <  Quadrangles_area < (reskm*(1+rdev_scale))**2

rzoom_fig = 5


if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    if len(argv) != 4:
        print('Usage: '+argv[0]+' <file_mojito.nc> <records to use (C), comma-separated> <basic_res_km>')
        exit(0)

    cf_nc_in = argv[1]
    lstrec = argv[2]
    creskm = argv[3]
    reskm = float(creskm)

    vcrec = split(',',lstrec)
    Nrec  = len(vcrec)
    vRec = np.array( [ int(vcrec[i]) for i in range(Nrec) ], dtype=int )
    del vcrec

    mjt.chck4f(cf_nc_in)

    if not path.exists('./npz'): mkdir('./npz')
    if iplot>0:
        cfdir = './figs/quadgener'
        if not path.exists(cfdir): mkdir(cfdir)

    rtolQuadA = 0.75 * reskm/20. ; # +- tolerance in [km] to accept a given scale. Ex: average scale of quadrangle = 15.9 km is accepted for 15 km !!

    print('\n *** Allowed deviation from '+creskm+' km for the mean scale of constructed quads (i.e. `sqrt(mean(Quad_areas))`) = ',rtolQuadA,'km')
    
    #########################################################################################################
    
    rQarea_min, rQarea_max = (reskm*(1.-rdev_scale))**2, (reskm*(1.+rdev_scale))**2
    print('\n *** Will retain quadrangles with an area comprised between '+str(round(rQarea_min,1))+' km^2 and '+str(round(rQarea_max,1))+' km^2\n')

    # Loading the data for the 2 selected records:
    Nt, nBmax, corigin, lTimePos = mjt.GetDimNCdataMJT( cf_nc_in )

    if lTimePos: print(' *** There is the "time_pos" data in input netCDF file! => gonna use it!')
    
    # We need a string to name the output files:
    if corigin == 'RGPS':
        cfstr = 'RGPS_'+split('_', path.basename(cf_nc_in))[2] ; # Basically the name of the batch
        
    elif split('_',corigin)[0] == 'NEMO-SI3':
        sl = split('_', path.basename(cf_nc_in))[0:5] ; # Basically the name of the batch
        cfstr = sl[0]+'_'+sl[1]+'_'+sl[2]+'_'+sl[4]
    else:
        print('FIXME for corigin = '+corigin+' !!!'); exit(0)
    
    if np.any(vRec>=Nt):
        print('ERROR: some of the specified records # are >= '+str(Nt)+'  !'); exit(0)
    #
    vdate = np.zeros( Nrec,  dtype=int )
    vIDs  = np.zeros( nBmax, dtype=int )
    xPosG = np.zeros( (Nrec,nBmax,2) )
    xPosC = np.zeros( (Nrec,nBmax,2) )
    pmsk  = np.zeros( (Nrec,nBmax), dtype='i1' )
    if lTimePos:
        timePos = np.zeros( (Nrec,nBmax), dtype=int )        
    #
    jr = 0
    for jrec in vRec[:]:
        if lTimePos:
            vdate[jr], zIDs, xPosG[jr,:,:], xPosC[jr,:,:], pmsk[jr,:], timePos[jr,:] = mjt.LoadNCdataMJT( cf_nc_in, krec=jrec, lmask=True, lGetTimePos=True  )
        else:
            vdate[jr], zIDs, xPosG[jr,:,:], xPosC[jr,:,:], pmsk[jr,:]                = mjt.LoadNCdataMJT( cf_nc_in, krec=jrec, lmask=True, lGetTimePos=False )
        #
        print( ' * jrec = ',jrec, ', mean date =',mjt.epoch2clock(vdate[jr]))    
        if jr==0:
            vIDs[:] = zIDs[:]
        else:
            if np.sum(zIDs[:]-vIDs[:])!=0:
                print('ERROR: ID fuck up in input file!') ; exit(0)
        jr=jr+1

    
    # Need some calendar info:
    NbDays = int( (vdate[1] - vdate[0]) / (3600.*24.) )
    cdt1 = mjt.epoch2clock(vdate[0] )
    cdt2 = mjt.epoch2clock(vdate[-1])

    print('    *  start and End dates => '+cdt1+' -- '+cdt2,' | number of buoys =>',np.sum(pmsk[0,:]), np.sum(pmsk[1,:]))
    print('        ==> nb of days =', NbDays)


    # STUPID: #fixme
    zXY   = np.zeros( (nBmax,2,Nrec) )
    zXY[:,0,:] = xPosC[:,:,1].T
    zXY[:,1,:] = xPosC[:,:,0].T
    zGC   = np.zeros( (nBmax,2,Nrec) )
    zGC[:,0,:] = xPosG[:,:,1].T
    zGC[:,1,:] = xPosG[:,:,0].T

    mask = np.zeros( nBmax      , dtype='i1') + 1  ; # Mask to for "deleted" points (to cancel)    
    ztim = np.zeros((nBmax,Nrec), dtype=int )
    zPnm = np.array( [ str(i) for i in vIDs ], dtype='U32' ) ; # Name for each point, based on 1st record...

    if lTimePos:
        ztim[:,:] = timePos[:,:].T
    else:
        for jp in range(nBmax): ztim[jp,:] = vdate[:]

    for jr in range(Nrec):
        print('\n   * Record #'+str(vRec[jr])+':')

        # Translate rji,rjj from TracIce to lon,lat and x,y:
        #zGC[:,:,jr], zXY[:,:,jr] = mjt.rJIrJJtoCoord( xJJs[:,jr], xJIs[:,jr], xIDs[:,jr], xlon_t, xlon_u, xlat_t, xlat_v )

        # Get rid of points to close to land (shrinks arrays!):
        mask[:] = mjt.MaskCoastal( zGC[:,:,jr], mask=mask[:], rMinDistFromLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc )

    # How many points left after elimination of buoys that get too close to land (at any record):
    NbP  = np.sum(mask)
    NbRM = nBmax-NbP
    print('\n *** '+str(NbP)+' / '+str(nBmax)+' points survived the dist2coast test => ', str(NbRM)+' points to delete!')

    if NbRM>0:        
        zPnm, vIDs, zGC, zXY, ztim = mjt.ShrinkArrays( mask, zPnm, vIDs, zGC, zXY, ztim )

        
    if idebug>0:
        for jr in range(Nrec):
            print('\n  DEBUG *** Record #'+str(vRec[jr])+':')
            for jc in range(NbP):
                print( ' * #'+str(jc)+' => Name: "'+zPnm[jc]+'": ID=',vIDs[jc],', X=',zXY[jc,0,jr],', Y=',zXY[jc,1,jr],
                       ', lon=',zGC[jc,0,jr],', lat=',zGC[jc,1,jr], ', time=',mjt.epoch2clock(ztim[jc,jr]) )
                if str(vIDs[jc])!=zPnm[jc]:
                    print(' Fuck Up!!!! => vIDs[jc], zPnm[jc] =',vIDs[jc], zPnm[jc] ); exit(0)
        print('')

    
    cdate0  = str.replace( mjt.epoch2clock(vdate[0], precision='D'), '-', '')

    for jr in range(Nrec):
        jrec = vRec[jr]

        print('\n\n *** QUAD-MESH GENERATION => record #'+str(jrec)+': record = '+str(jr))

        cdats  = mjt.epoch2clock(vdate[jr])
        cdate  = str.replace( mjt.epoch2clock(vdate[jr], precision='D'), '-', '')

        cfbase = cfstr+'_'+cdate0+'t0_'+cdate+'_'+creskm+'km'
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
                # fed into Delaunay (zXY), if it is the case we have to shrink `zXY`, `vIDs`, `zPnm`
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
                    zPnm, vIDs, zGC, zXY, ztim = mjt.ShrinkArrays( mask, zPnm, vIDs, zGC, zXY, ztim )
                    del zPntIdx0, mask, ind2keep
                    NbP = NbPnew
                    print('      => done! Only '+str(NbP)+' points left in the game...')

                del PntIdxInUse, NbPnew

            ### while not lOK
            print('\n *** We have '+str(NbT)+' triangles!')

            
            # Conversion to the `Triangle` class:
            TRIAS = mjt.Triangle( zXY[:,:,jr], xTpnts, xNeighborIDs, vIDs, ztim[:,jr], zPnm, origin=corigin )
            del xTpnts, xNeighborIDs, TRI

            # Info on the triangles:
            if idebug>0:
                zlngth = np.mean( TRIAS.lengths(), axis=1)
                zml = np.mean(zlngth)
                print('   => mean length and stdev of the sides of all triangles =',zml,mjt.StdDev( zml, zlngth ), 'km')
                del zlngth, zml
            
            # Merge triangles into quadrangles:
            xQcoor, vPids, vTime, xQpnts, vQnam = mjt.Tri2Quad( TRIAS, anglRtri=(rTang_min,rTang_max),
                                                                ratioD=rdRatio_max, anglR=(rQang_min,rQang_max),
                                                                areaR=(rQarea_min,rQarea_max), idbglev=idebug )
            if len(xQpnts)<=0: exit(0)

            (NbQ,_) = np.shape(xQpnts)
            print('\n *** We have '+str(NbQ)+' quadrangles!')

            # Save the triangular mesh info:
            mjt.SaveClassPolygon( cf_npzT, TRIAS, ctype='T', origin=corigin )

            # To be used for other record, indices of Points to keep for Quads:
            _,ind2keep,_ = np.intersect1d(vIDs, vPids, return_indices=True); # retain only indices of `vIDs` that exist in `vPids`


            # Plots specific to `jrec==0`:
            if iplot>0:
                # We need to find a descent X and Y range for the figures:
                vrngX = mjt.roundAxisRange( TRIAS.PointXY[:,0], rndKM=50. )
                vrngY = mjt.roundAxisRange( TRIAS.PointXY[:,1], rndKM=50. )
                #
                if iplot>3:
                    # Show all initial points (out of TrackIce):
                    print('\n *** Launching initial cloud point plot!')
                    kk = mjt.ShowTQMesh( zXY_dbg[:,0], zXY_dbg[:,1], cfig=cfdir+'/fig01a_cloud_points_'+cfbase+'.png',
                                         ppntIDs=vPids_dbg, lGeoCoor=False, zoom=rzoom_fig,
                                         rangeX=vrngX, rangeY=vrngY )
                if iplot>2:
                    # Show triangles on a map:
                    print('\n *** Launching Triangle plot!')
                    kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig=cfdir+'/fig01_Mesh_Triangles_'+cfbase+'.png',
                                         ppntIDs=TRIAS.PointIDs,
                                         TriMesh=TRIAS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig,
                                         rangeX=vrngX, rangeY=vrngY )
                    
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS0 = mjt.Quadrangle( xQcoor, xQpnts, vPids, vTime, vQnam, date=cdats, origin=corigin )


            # Some info about the spatial scales of quadrangles:
            print('\n *** About our quadrangles:')
            zsides = QUADS0.lengths()
            zareas = QUADS0.area()
            rl_average_side = np.mean(zsides)
            rl_average_scal = np.mean( np.sqrt(zareas) )
            rl_average_area = np.mean(zareas) ; rl_stdev_area = mjt.StdDev(rl_average_area, zareas)
            print('    ==> average scale (sqrt[A]) is '+str(round(rl_average_scal,3))+' km')
            print('    ==> average side length is '+str(round(rl_average_side,3))+' km')
            print('    ==> average area is '+str(round(rl_average_area,1))+' km^2, StDev =',str(round(rl_stdev_area,1))+' km^2')
            del zareas, zsides
            if abs(rl_average_scal-reskm) > rtolQuadA:
                print(' ERROR: this avergae scale is too different from the '+creskm+' km expected!!!')
                exit(0)
            
            # Save the quadrangular mesh info:
            mjt.SaveClassPolygon( cf_npzQ, QUADS0, ctype='Q', origin=corigin )


            #######################################################################################################
        else:
            
            #  => jr >= 1 !!!

            print('\n *** Quad recycling for jr=',jr)

            # Recycling Quads found at 1st record (QUADS0): lilo
            xQcoor, vTime, xQpnts, vPids, vQnam, vQIDs  = mjt.RecycleQuads( zXY[:,:,jr], ztim[:,jr], vIDs, QUADS0,  iverbose=idebug )
            
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS = mjt.Quadrangle( xQcoor, xQpnts, vPids, vTime, vQnam, vQIDs=vQIDs, date=cdats, origin=corigin )

            # Save the quadrangular mesh info:
            mjt.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q', origin=corigin )

            del QUADS

        ### if jr == 0


        if iplot>0:

            TRI = mjt.LoadClassPolygon( cf_npzT, ctype='T' )
            QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )


            if iplot>1:
                # Show triangles together with the quadrangles on a map:
                print('\n *** Launching Triangle+Quad plot!')
                kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig=cfdir+'/fig02_Mesh_Quadrangles_'+cfbase+'.png',
                                     ppntIDs=TRI.PointIDs, TriMesh=TRI.MeshVrtcPntIdx,
                                     pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshVrtcPntIdx,
                                     lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
            if iplot>2:
                # Show only points composing the quadrangles:
                kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig=cfdir+'/fig03a_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                     ppntIDs=QUA.PointIDs, lGeoCoor=False, zoom=rzoom_fig )
            
            # Show only the quads with only the points that define them:
            print('\n *** Launching Quad-only plot!')
            kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig=cfdir+'/fig03_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                 ppntIDs=QUA.PointIDs, QuadMesh=QUA.MeshVrtcPntIdx, qIDs=QUA.QuadIDs,
                                 lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )

    ### for jr in range(Nrec):
    print('\n --- Over!\n')

