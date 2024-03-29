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
iplot  = 2 ; # Create figures to see what we are doing...

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

rtolQuadA        = 0.5 ; # +- tolerance in [km] to accept a given scale. Ex: 15.19 km is accepted for 15 km !!


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

    lNeed4SubSamp = ( reskm > 12. )
    
    vcrec = split(',',lstrec)
    Nrec  = len(vcrec)
    vRec = np.array( [ int(vcrec[i]) for i in range(Nrec) ], dtype=int )
    del vcrec

    mjt.chck4f(cf_nc_in)

    if not path.exists('./npz'): mkdir('./npz')
    if iplot>0:
        cfdir = './figs/quadgener'
        if not path.exists(cfdir): mkdir(cfdir)

    #########################################################################################################


    if reskm>7.5 and reskm<12.5:
        rQarea_min =  60. ; # min area allowed for Quadrangle [km^2]
        rQarea_max = 200. ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:1
        #
    elif reskm>17.5 and reskm<22.5:
        rQarea_min = 300.
        rQarea_max = 600.
        #
    elif reskm>20 and reskm<30:
        rQarea_min = 350.
        rQarea_max = 950.
        #
    elif reskm>30 and reskm<40:
        rQarea_min =  800.
        rQarea_max = 1800.
        #
    elif reskm>50 and reskm<70:
        rQarea_min = 2000.
        rQarea_max = 5000.        
        #
    elif reskm>200 and reskm<400:
        rQarea_min =  10000.
        rQarea_max = 50000.       
        #
    else:
        print('ERROR: we do not know what to do with resolution `reskm` =', reskm) ; exit(0)
        #
        #rQarea_max = 18000. ; rzoom_fig = 2  ; # max area allowed for Quadrangle [km^2] VALID for NANUK4 HSS:10



    print(' *** Specified spatial scale is '+creskm+' km')
    if lNeed4SubSamp:
        print('   => the cloud of point will be sub-sampled!')
    print('   => will retain quadrangles with an area: ',rQarea_min,'km^2 < A <',rQarea_max,'km^2')
        
        

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
        mask[:] = mjt.MaskCoastal( zGC[:,:,jr], mask=mask[:], rMinDistLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc )

    # How many points left after elimination of buoys that get too close to land (at any record):
    NbP1 = np.sum(mask)
    print('\n *** '+str(NbP1)+' / '+str(nBmax)+' points survived the dist2coast test => ', str(nBmax-NbP1)+' points to delete!')
    
    NbP = NbP1

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
        cf_npzT = './npz/T-mesh_'+cfbase+'.npz'

        
        if jr == 0:

            #lili:
            l_happy = False
            itt     = 0
            rfcorr  = 0.7
            rdev    = 1.
            while not l_happy:
                itt = itt + 1

                if lNeed4SubSamp:
                    
                    # SUB-SAMPLING
                    #  From experience we have to pick a scale significantly smaller that that of the desired
                    #  quadrangles in order to obtain quadrangle of the correct size!

                    if reskm==20.:
                        rd_ss = 15.1
                    else:
                        rd_ss = rfcorr * reskm ; # Correct `reskm` to get closer to requested radius (based on QUADs to be generated)
    
                    print('\n *** Applying spatial sub-sampling with radius: '+str(round(rd_ss,2))+'km')
                    #print('LOLO: shape(zXY) =',np.shape(zXY))
                    NbPss, zXYss, idxKeep = mjt.SubSampCloud( rd_ss, zXY[:,:,jr] )
                    #exit(0)
                else:
                    # No sub-sampling:
                    NbPss = NbP                    
                    zXYss = zXY[:,:,jr].copy()
                    idxKeep = np.arange(NbP)
                    l_happy = True

                
                print('\n *** Delaunay triangulation for 1st record! pass #',itt)                    
    
                # Generating triangular meshes out of the cloud of points for 1st record:
                lOK = False
                while not lOK:
    
                    TRI = Delaunay( zXYss )
                    
                    xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*
                    (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles
                    xNeighborIDs = TRI.neighbors.copy() ;  # shape = (Nbt,3)
    
                    # For some reasons, here, `xTpnts` can involve less points than the number of points
                    # fed into Delaunay (zXYss), if it is the case we have to shrink `zXY`, `vIDs`, `zPnm`
                    # accordingly...
                    PntIdxInUse = np.unique( xTpnts.flatten() )
                    NbPnew      = len(PntIdxInUse)
                    lOK = ( NbPnew == NbPss )
    
                    if not lOK:
                        if NbPnew > NbPss:
                            print('ERROR: `NbPnew > NbPss` !!!'); exit(0)
                        print('\n *** Need to cancel the '+str(NbPss-NbPnew)+' points that vanished in Delaunay triangulation!')
                        mask  = np.zeros(  NbPss, dtype=int )
                        zPntIdx0 = np.arange( NbPss, dtype=int )
                        _,ind2keep,_ = np.intersect1d(zPntIdx0, PntIdxInUse, return_indices=True); # retain only indices of `zPntIdx0` that exist in `PntIdxInUse`
                        mask[ind2keep] = 1
                        if np.sum(mask) != NbPnew:
                            print('ERROR: `np.sum(mask) != NbPnew` !!!'); exit(0)
                        zPnm, vIDs, zGC, zXY, ztim = mjt.ShrinkArrays( mask, zPnm, vIDs, zGC, zXY, ztim )
                        del zPntIdx0, mask, ind2keep
                        NbPss = NbPnew
                        print('      => done! Only '+str(NbPss)+' points left in the game...')
    
                    del PntIdxInUse, NbPnew
    
                ### while not lOK
                print('\n *** We have '+str(NbT)+' triangles! pass #',itt)
    
                
                # Conversion to the `Triangle` class:
                TRIAS = mjt.Triangle( zXYss, xTpnts, xNeighborIDs, vIDs[idxKeep], ztim[idxKeep,jr], zPnm[idxKeep], origin=corigin )
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
                #if len(xQpnts)<=0: exit(0)
                #(NbQ,_) = np.shape(xQpnts)
                #print('\n *** We have '+str(NbQ)+' quadrangles!')
                l_someQuads = (len(xQpnts)>0)
                #if len(xQpnts)<=0: exit(0)
                
                if l_someQuads:
                    (NbQ,_) = np.shape(xQpnts)
                    print('\n *** We have '+str(NbQ)+' quadrangles! pass #',itt)

                    # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
                    QUADS0 = mjt.Quadrangle( xQcoor, xQpnts, vPids, vTime, vQnam, date=cdats, origin=corigin )
                        
                    zsides = QUADS0.lengths()
                    zareas = QUADS0.area()
                    rl_average_side = np.mean(zsides)
                    rl_average_scal = np.mean( np.sqrt(zareas) )
                    rl_average_area = np.mean(zareas) ; rl_stdev_area = mjt.StdDev(rl_average_area, zareas)
                    print('    ==> average scale (sqrt[A]) is '+str(round(rl_average_scal,3))+' km')
                    print('    ==> average side length is '+str(round(rl_average_side,3))+' km')
                    print('    ==> average area is '+str(round(rl_average_area,1))+' km^2, StDev =',str(round(rl_stdev_area,1))+' km^2')
                    del zareas, zsides

                    if lNeed4SubSamp:
                        rdev_old = rdev
                        rdev     = rl_average_scal - reskm
                        l_happy  = ( abs(rdev) < rtolQuadA ) ; # average quadrangle scale (root square of area) is close to expected nominal scale

                        # Surrender block:
                        #if not l_happy:
                        #    # If we are at the nominal scale:
                        #    l_happy = ( int(reskm) == int(rd_nom_data) and itt>4 )
                        #if not l_happy:
                        #    l_happy = (itt>7)

                        # Refining:
                        if not l_happy:                                
                            # Linear fit of actual correction as a function of `reskm`
                            #rfc = 0.008*reskm + 0.56
                            #rfc = 0.008*reskm + 0.7
                            rfc = 0.008*reskm + 0.6                
                            if itt==1: ralpha = (1.-rfc) / rdev ; # equivalent to a correction of `rfc`
                            if itt>1 and copysign(1,rdev) == -copysign(1,rdev_old):
                                ralpha = ralpha/1.3 ; # change of sign of deviation => we decrease alpha!
                            # will go for a next round with a correction factor becoming increasingly smaller than 1:
                            rfcorr = min( max(0.6 , rfcorr - ralpha * rdev ) , 0.95 )
                            #rfcorr = min( max(0.5 , rfcorr - ralpha * rdev ) , 0.95 )
                            print(' +++++ NEW itteration with: ralpha, rfcorr = ',ralpha, rfcorr)
                        #
                else:
                    l_happy = True ; # we coud not build any QUAD so we move on anyways...
                #
                ### if l_someQuads
                    
                    

            ### while not l_happy

            if iplot>0:
                # We need to find a descent X and Y range for the figures:
                vrngX = mjt.roundAxisRange( TRIAS.PointXY[:,0], rndKM=50. )
                vrngY = mjt.roundAxisRange( TRIAS.PointXY[:,1], rndKM=50. )

            
            if iplot >1 and lNeed4SubSamp:
                # Shows the cloud of buoys (with buoys' IDs) on the Cartesian plane (km)
                # After and before subsampling
                kk = mjt.ShowTQMesh( zXY[:,0,0], zXY[:,1,0], cfig=cfdir+'/00_'+cfbase+'_Original.png',
                                     ppntIDs=vIDs[:], lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
                kk = mjt.ShowTQMesh( zXYss[:,0], zXYss[:,1], cfig=cfdir+'/00_'+cfbase+'_SubSamp.png',
                                     ppntIDs=vIDs[idxKeep], lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )

            
            # Save the triangular mesh info:
            mjt.SaveClassPolygon( cf_npzT, TRIAS, ctype='T', origin=corigin )
            
            # To be used for other record, indices of Points to keep for Quads:
            _,ind2keep,_ = np.intersect1d(vIDs, vPids, return_indices=True); # retain only indices of `vIDs` that exist in `vPids`


            # Plots specific to `jrec==0`:
            if iplot>2:
                # Show all initial points (out of TrackIce):
                print('\n *** Launching initial cloud point plot!')
                kk = mjt.ShowTQMesh( zXY_dbg[:,0], zXY_dbg[:,1], cfig=cfdir+'/fig01a_cloud_points_'+cfbase+'.png',
                                     ppntIDs=vPids_dbg, lGeoCoor=False, zoom=rzoom_fig,
                                     rangeX=vrngX, rangeY=vrngY )
            if iplot>1:
                # Show triangles on a map:
                print('\n *** Launching Triangle plot!')
                kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig=cfdir+'/fig01_Mesh_Triangles_'+cfbase+'.png',
                                     ppntIDs=TRIAS.PointIDs,
                                     TriMesh=TRIAS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig,
                                     rangeX=vrngX, rangeY=vrngY )
                    
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS0 = mjt.Quadrangle( xQcoor, xQpnts, vPids, vTime, vQnam, date=cdats, origin=corigin )

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

            QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )

            if iplot>1 and path.exists(cf_npzT):
                TRI = mjt.LoadClassPolygon( cf_npzT, ctype='T' )
                # Show triangles together with the quadrangles on a map:
                print('\n *** Launching Triangle+Quad plot!')
                kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig=cfdir+'/fig02_Mesh_Quadrangles_'+cfbase+'.png',
                                     ppntIDs=TRI.PointIDs, TriMesh=TRI.MeshVrtcPntIdx,
                                     pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshVrtcPntIdx,
                                     lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
            if iplot>2:
                # Show only points composing the quadrangles:
                kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig=cfdir+'/fig03a_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                     ppntIDs=QUA.PointIDs, lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
            
            # Show only the quads with only the points that define them:
            print('\n *** Launching Quad-only plot!')
            kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig=cfdir+'/fig03_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                 ppntIDs=QUA.PointIDs, QuadMesh=QUA.MeshVrtcPntIdx, qIDs=QUA.QuadIDs,
                                 lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )

    ### for jr in range(Nrec):
    print('\n --- Over!\n')

