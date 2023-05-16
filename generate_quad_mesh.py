#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#  INPUT DATA: a `nc` file created with `ncSaveCloudBuoys()` of `mojito/ncio.py` !
#
#    L. Brodeau, 2023
#
#
#
# TO DO:
#
##################################################################

from sys import argv, exit
from os import path, environ, mkdir, makedirs
import numpy as np
from re import split
from scipy.spatial import Delaunay

import mojito   as mjt
from mojito import config as cfg

idebug = 0

iplot  = 0 ; # Create figures to see what we are doing...

rzoom_fig = 5

lapplyDistCoast = False

lExportCloudPoints = False ; # in case we want to save the location of valid/selected quads in a netCDF file with `ncSaveCloudBuoys`
#                           # => in order to seed from it!

if __name__ == '__main__':    
    
    crd_ss = None





    
    if not len(argv) in [5,6]:
        print('Usage: '+argv[0]+' <file_mojito.nc> <records to use (C), comma-separated> <basic_res_km> <mode (rgps,model,xlose)> (<rd_ss_km>)')
        exit(0)






        
    cf_nc_in = argv[1]
    lstrec = argv[2]
    creskm = argv[3]
    reskm = int(creskm)
    quality_mode = argv[4]

    if not quality_mode  in ['thorough','model','rgps','xlose']:
        print('ERROR => unknow mode: '+mode+'!') ; exit(0)
    
    if len(argv)==6:
        crd_ss = argv[5]

    
    vcrec = split(',',lstrec)
    Nrec  = len(vcrec)
    vRec = np.array( [ int(vcrec[i]) for i in range(Nrec) ], dtype=int )
    del vcrec

    mjt.chck4f(cf_nc_in)
    
    if not path.exists('./npz'): mkdir('./npz')
    if iplot>0:
        fdir = './figs/quadgener/'+str(reskm)+'km'
        makedirs( fdir, exist_ok=True )
    
    # Loading the data for the 2 selected records:
    Nt, NbP, corigin, lTimePos = mjt.GetDimNCdataMJT( cf_nc_in )

    if lTimePos: print(' *** There is the "time_pos" data in input netCDF file! => gonna use it!')

    if not corigin in ['RGPS', 'model']
        print('ERROR: data origin "'+corigin+'" is unknown !!! => Fix me!!!')
        exit(0)
        
    kl = cfg.initialize( mode=quality_mode )
    kk = cfg.updateConfig4Scale( reskm, mode=quality_mode )

    zA_ideal = reskm*reskm ; # ideal Area of quadrangles to build!
    
    print('\n *** Allowed deviation from '+creskm+' km for the MEAN scale of constructed quads (i.e. `sqrt(mean(Quad_areas))`) = ',cfg.rc_tolQuadA,'km')
    print(' *** Max time deviation accepted for vertices: `rc_t_dev_cancel` =',cfg.rc_t_dev_cancel,'s')

    print('\n *** Will only retain quadrangles with an area comprised between '
          +str(round(cfg.rc_Qarea_min,1))+' km^2 and '+str(round(cfg.rc_Qarea_max,1))+' km^2, Nominal='+str(zA_ideal)+'km^2 \n')

    
    # We need a string to name the output files:
    kdt = -4
    #if creskm=='10': kdt = -3
    #
    if corigin == 'RGPS':
        cfstr = 'RGPS_'+split('_', path.basename(cf_nc_in))[2] ; # Basically the name of the batch
        cdtbin = '_'+split('_', path.basename(cf_nc_in))[kdt] ; print('         > cdtbin =',cdtbin)
        #
    elif split('_',corigin)[0] == 'NEMO-SI3':
        sl = split('_', path.basename(cf_nc_in))[0:5] ; # Basically the name of the batch
        cfstr = sl[0]+'_'+sl[1]+'_'+sl[2]+'_'+sl[4]
        cdtbin = '_'+split('_', path.basename(cf_nc_in))[kdt] ; print('         > cdtbin =',cdtbin)
    if not cdtbin[1:3] in ['dt']:
        print('WARNING: we could not figure out `cdtbin`!')
        cdtbin = '_NoBin'
    cfstr += cdtbin
    

    if np.any(vRec>=Nt):
        print('ERROR: some of the specified records # are >= '+str(Nt)+'  !'); exit(0)
    #
    vdate = np.zeros( Nrec,  dtype=int )
    vIDs  = np.zeros( NbP, dtype=int )
    zGC   = np.zeros( (NbP,2,Nrec) )
    zXY   = np.zeros( (NbP,2,Nrec) )
    zmsk  = np.zeros( (NbP,Nrec), dtype='i1' )
    ztim  = np.zeros( (NbP,Nrec), dtype=int )

    jr = 0
    for jrec in vRec[:]:
        if lTimePos:
            vdate[jr], zIDs, zGC[:,:,jr], zXY[:,:,jr], zmsk[:,jr], ztim[:,jr] = mjt.LoadNCdataMJT( cf_nc_in, krec=jrec, lmask=True,
                                                                                                          lGetTimePos=True, convention='F' )
        else:
            vdate[jr], zIDs, zGC[:,:,jr], zXY[:,:,jr], zmsk[:,jr]             = mjt.LoadNCdataMJT( cf_nc_in, krec=jrec, lmask=True,
                                                                                                          lGetTimePos=False, convention='F' )
        #
        
        print( ' * jrec = ',jrec, ', mean date =',mjt.epoch2clock(vdate[jr]))    
        if jr==0:
            vIDs[:] = zIDs[:]
        else:
            if np.sum(zIDs[:]-vIDs[:])!=0:
                print('ERROR: ID fuck up in input file!') ; exit(0)
        jr=jr+1

        
    # Need some calendar info:
    NbDays = int( (vdate[1] - vdate[0]) / cfg.rc_day2sec )
    cdt1 = mjt.epoch2clock(vdate[0] )
    cdt2 = mjt.epoch2clock(vdate[-1])

    print('    *  start and End dates => '+cdt1+' -- '+cdt2,' | number of buoys =>',np.sum(zmsk[:,0]), np.sum(zmsk[:,1]))
    print('        ==> nb of days =', NbDays)

    # Some sanity controls on time
    if Nt==2:
        zdt = vdate[1] - vdate[0]
        print('     => based on these 2 dates the "mean" `dt` =',round(zdt/(3600*24),1),' days')
        # ===> dt_Nmnl, max_dev_dt_Nmnl

    
    if lTimePos and Nt==2:
        # Time control:
        zzt = ztim[:,1] - ztim[:,0]
        if np.any( zzt == 0. ):
            print('ERROR: some identical times in the 1st and 2nd records!!!')
            (idxFU,) = np.where( zzt == 0. )
            print('  => for '+str(len(idxFU))+' points!')
            exit(0)
        del zzt
    
    mask = np.zeros( NbP      , dtype='i1') + 1  ; # Mask to for "deleted" points (to cancel)    
    zPnm = np.array( [ str(i) for i in vIDs ], dtype='U32' ) ; # Name for each point, based on 1st record...

    if not lTimePos:
        for jp in range(NbP): ztim[jp,:] = vdate[:]

    if lapplyDistCoast:
        for jr in range(Nrec):
            print('\n   * Record #'+str(vRec[jr])+':')
    
            # Translate rji,rjj from TracIce to lon,lat and x,y:
            #zGC[:,:,jr], zXY[:,:,jr] = mjt.rJIrJJtoCoord( xJJs[:,jr], xJIs[:,jr], xIDs[:,jr], xlon_t, xlon_u, xlat_t, xlat_v )
    
            # Get rid of points to close to land (shrinks arrays!):
            mask[:] = mjt.MaskCoastal( zGC[:,:,jr], mask=mask[:], rMinDistLand=cfg.nc_MinDistFromLand, fNCdist2coast=cfg.fdist2coast_nc )
    
        # How many points left after elimination of buoys that get too close to land (at any record):
        NbPn = np.sum(mask)
        NbRM = NbP-NbPn
        print('\n *** '+str(NbPn)+' / '+str(NbP)+' points survived the dist2coast test => ', str(NbRM)+' points to delete!')
    
        if NbRM>0:
            zPnm, vIDs, zGC, zXY, ztim = mjt.ShrinkArrays( mask, zPnm, vIDs, zGC, zXY, ztim )
            NbP = NbPn

        
    if idebug>0:
        for jr in range(Nrec):
            print('\n  DEBUG *** Record #'+str(vRec[jr])+':')
            for jc in range(NbP):
                print( ' * #'+str(jc)+' => Name: "'+zPnm[jc]+'": ID=',vIDs[jc],', X=',zXY[jc,0,jr],', Y=',zXY[jc,1,jr],
                       ', lon=',zGC[jc,0,jr],', lat=',zGC[jc,1,jr], ', time=',mjt.epoch2clock(ztim[jc,jr]) )
                if str(vIDs[jc])!=zPnm[jc]:
                    print(' Fuck Up!!!! => vIDs[jc], zPnm[jc] =',vIDs[jc], zPnm[jc] ); exit(0)
        print('')

    
    cdate0  = str.replace( mjt.epoch2clock(vdate[0], precision='m'), '-', '')
    cdate0   = str.replace( cdate0, '_', '-')
    cdate0   = str.replace( cdate0, ':', 'h')


    if NbP<4:
        print('  \n *** Only '+str(NbP)+' points in the cloud! => exiting!!!');  exit(0)    

    
    for jr in range(Nrec):
        jrec = vRec[jr]

        print('\n\n *** QUAD-MESH GENERATION => record #'+str(jrec)+': record = '+str(jr))

        cdats  = mjt.epoch2clock(vdate[jr])
        cdate  = str.replace( mjt.epoch2clock(vdate[jr], precision='m'), '-', '')
        cdate   = str.replace( cdate, '_', '-')
        cdate   = str.replace( cdate, ':', 'h')

        cfbase = cfstr+'_'+cdate0+'t0_'+cdate+'_'+creskm+'km'
        if crd_ss:
            cfbase = cfstr+'_'+cdate0+'t0_'+cdate+'_'+crd_ss+'-'+creskm+'km'
        
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
            xPxy, vPids, vTime, xQpnts, vQnam = mjt.Tri2Quad( TRIAS, anglRtri=(cfg.rc_Tang_min,cfg.rc_Tang_max),
                                                              ratioD=cfg.rc_dRatio_max, anglR=(cfg.rc_Qang_min,cfg.rc_Qang_max),
                                                              areaR=(cfg.rc_Qarea_min,cfg.rc_Qarea_max), areaIdeal=zA_ideal,
                                                              idbglev=idebug )

            lquit = ( np.shape(xPxy)==(0,) or np.shape(xQpnts)==(0,) )
            if not lquit:
                (nbP,_), (nbQ,_) = np.shape(xPxy), np.shape(xQpnts)                
                lquit = (nbQ<=0)
            if lquit:
                print('  \n *** 0 quads from triangles!!! => exiting!!!');  exit(0)
                
            print('\n *** After `Tri2Quad()`: we have '+str(nbQ)+' quadrangles out of '+str(nbP)+' points!')


            if corigin == 'RGPS':
                # Check if time deviation accros the 4 vertices of a quadrangle is too large, and cancel it if so...
                xPxy, vPids, vTime, xQpnts, vQnam = mjt.QuadVrtcTDev( xPxy, vPids, vTime, xQpnts, vQnam, cfg.rc_t_dev_cancel )
                (nbP,_), (nbQ,_) = np.shape(xPxy), np.shape(xQpnts)
                print('\n *** After vertices time check: we have '+str(nbQ)+' quadrangles out of '+str(nbP)+' points')

                
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
                if iplot>2:
                    # Show triangles on a map:
                    print('\n *** Launching Triangle plot!')
                    kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig=fdir+'/fig01_Mesh_Triangles_'+cfbase+'.png',
                                         ppntIDs=TRIAS.PointIDs,
                                         TriMesh=TRIAS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig,
                                         rangeX=vrngX, rangeY=vrngY )
                    
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS0 = mjt.Quadrangle( xPxy, xQpnts, vPids, vTime, vQnam, date=cdats, origin=corigin, reskm_nmnl=reskm )

            #exit(0); #lili
            # Some info about the spatial scales of quadrangles:
            print('\n *** About our quadrangles:')
            zsides = QUADS0.lengths()
            zareas = QUADS0.area()
            zscale = np.sqrt(zareas)
            rl_average_side = np.mean(zsides)
            rl_average_scal = np.mean(zscale)
            rl_stdev_scal   = mjt.StdDev(rl_average_scal, zscale)
            rl_average_area = np.mean(zareas)
            print('    ==> scale: mean, StDev, min, max =',round(rl_average_scal,3), round(rl_stdev_scal,2), round(np.min(zscale),2), round(np.max(zscale),3),' km')
            print('    ==> average side length is '+str(round(rl_average_side,3))+' km')
            print('    ==> average area is '+str(round(rl_average_area,1))+' km^2')
            del zareas, zsides, zscale
            zdev = abs(rl_average_scal-reskm)
            if zdev > cfg.rc_tolQuadA:
                print(' ERROR: the mean scale is too different from the '+creskm
                      +'km expected!!! (dev.=',round(zdev,2),' tol. = '+str(cfg.rc_tolQuadA)+'km)')
                exit(0)
            
            # Save the quadrangular mesh info:
            mjt.SaveClassPolygon( cf_npzQ, QUADS0, ctype='Q', origin=corigin, reskm_nmnl=reskm )


            #
            if lExportCloudPoints:  
                # To be save in a NC file:
                zpid  = np.zeros(nbP, dtype=int)
                zPxy  = np.zeros((nbP,2,Nrec))
                zPgc  = np.zeros((nbP,2,Nrec))
                zmask = np.zeros((nbP,Nrec))
                ztime = np.zeros((nbP,Nrec))
                zPxy[:,:,jr] = xPxy[:,:]
                for jb in range(nbP):
                    [zx,zy] = xPxy[jb,:]
                    #print(' zx,zy =',zx,zy)
                    (idx,) = np.where( (np.round(zXY[:,0,jr],6) == round(zx,6)) & (np.round(zXY[:,1,jr],6) == round(zy,6)) )
                    #print('idx =',idx)
                    if len(idx)!=1:
                        print('ERROR: `len(idx)!=1` !'); exit(0)
                    zPgc[jb,:,jr] =  zGC[idx[0],:,jr]
                    zpid[jb]      = vIDs[idx[0]]
                    zmask[jb,jr] = zmsk[idx[0],jr]
                    ztime[jb,jr] =  ztim[idx[0],jr]
                    #print('zXY[idx,:,jr], xPxy[jb,:]',zXY[idx,:,jr], xPxy[jb,:])
                    
            #exit(0)
                    

                
            
            
            #######################################################################################################
        else:
            
            #  => jr >= 1 !!!

            print('\n *** Quad recycling for jr=',jr)

            # Recycling Quads found at 1st record (QUADS0): lilo
            xPxy, vTime, xQpnts, vPids, vQnam, vQIDs  = mjt.RecycleQuads( zXY[:,:,jr], ztim[:,jr], vIDs, QUADS0,  iverbose=idebug )
            
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS = mjt.Quadrangle( xPxy, xQpnts, vPids, vTime, vQnam, vQIDs=vQIDs, date=cdats, origin=corigin, reskm_nmnl=reskm )

            # Save the quadrangular mesh info:
            mjt.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q', origin=corigin, reskm_nmnl=reskm )
            
            del QUADS

            
            if lExportCloudPoints:
                # To be save in a NC file:
                zPxy[:,:,jr] = xPxy[:,:]
                for jb in range(nbP):
                    [zx,zy] = xPxy[jb,:]
                    (idx,) = np.where( (np.round(zXY[:,0,jr],6) == round(zx,6)) & (np.round(zXY[:,1,jr],6) == round(zy,6)) )
                    if len(idx)!=1:
                        print('ERROR: `len(idx)!=1` !'); exit(0)
                    zPgc[jb,:,jr] =  zGC[idx[0],:,jr]
                    zpid[jb]      = vIDs[idx[0]]
                    zmask[jb,jr]  = zmsk[idx[0],jr]
                    ztime[jb,jr]  =  ztim[idx[0],jr]





            
        ### if jr == 0



        


        
        if iplot>0:

            TRI = mjt.LoadClassPolygon( cf_npzT, ctype='T' )
            QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )


            if iplot>1:
                # Show triangles together with the quadrangles on a map:
                print('\n *** Launching Triangle+Quad plot!')
                kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig=fdir+'/fig02_Mesh_Quadrangles_'+cfbase+'.png',
                                     ppntIDs=TRI.PointIDs, TriMesh=TRI.MeshVrtcPntIdx,
                                     pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshVrtcPntIdx,
                                     lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
            if iplot>2:
                # Show only points composing the quadrangles:
                kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig=fdir+'/fig03a_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                     ppntIDs=QUA.PointIDs, lGeoCoor=False, zoom=rzoom_fig )
            
            # Show only the quads with only the points that define them:
            print('\n *** Launching Quad-only plot!')
            kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig=fdir+'/fig03_Mesh_Points4Quadrangles_'+cfbase+'.png',
                                 ppntIDs=QUA.PointIDs, QuadMesh=QUA.MeshVrtcPntIdx, qIDs=QUA.QuadIDs,
                                 lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )


            #kk = mjt.ShowQuads( QUA.PointXY[QUA.MeshVrtcPntIdx,:], cfig=fdir+'/fig04_Quadrangles_'+cfbase+'.png', rangeX=vrngX, rangeY=vrngY )
            
    ### for jr in range(Nrec):





    if lExportCloudPoints:
        cf_nc_out = str.replace( cf_nc_in, '.nc', '_postQG.nc' )

        zYkm = np.zeros( (Nrec,nbP) )
        zXkm = np.zeros( (Nrec,nbP) )
        zlat = np.zeros( (Nrec,nbP) )
        zlon = np.zeros( (Nrec,nbP) )

        zYkm[:,:] = zPxy[:,1,:].T
        zXkm[:,:] = zPxy[:,0,:].T
        zlat[:,:] = zPgc[:,1,:].T
        zlon[:,:] = zPgc[:,0,:].T

        print('   * saving '+cf_nc_out)
        
        kk = mjt.ncSaveCloudBuoys( cf_nc_out, vdate, zpid, zYkm, zXkm, zlat, zlon, mask=zmask.T,
                                   xtime=ztime.T, fillVal=mjt.FillValue, corigin=corigin )

        if iplot>0:
            # Ploting what we saved in NC file:
            for jr in range(Nrec):
                cfig = str.replace( cf_nc_out, '.nc', '_postQG'+'_rec%3.3i'%(jr)+'.png' )
                cfig = str.replace( cfig , './nc', fdir )
                kk = mjt.ShowBuoysMap( vdate[jr], zlon[jr,:], zlat[jr,:], pvIDs=vIDs[:], cfig=cfig,
                                       cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )


        
    
    print('\n --- Over!\n')




    
