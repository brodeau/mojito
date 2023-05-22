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
from mojito.util import epoch2clock as e2c

idebug = 0

iplot = 0 ; # Create figures to see what we are doing...

rzoom_fig = 5

lExportCloudPoints = False ; # in case we want to save the location of valid/selected quads in a netCDF file with `ncSaveCloudBuoys`
#                           # => in order to seed from it!

mode_ctl_vrtc = 'maxdiff'


def _formDate_( idate ):
    cdate = str.replace( e2c(idate, precision='m'), '-', '')
    cdate = str.replace( cdate, '_', '-')
    cdate = str.replace( cdate, ':', 'h')
    return cdate


def _init4rec_( kr, pRec, pdate, cdt0, cfs, crkm, css ):
    krec = pRec[kr]
    print('\n\n *** QUAD-MESH GENERATION => record #'+str(krec)+': rec_itt = '+str(kr))
    #
    cdts = e2c(pdate[kr])
    cdte = _formDate_( pdate[kr] )
    cfbs = cfs+'_'+cdt0+'t0_'+cdte+'_'+crkm+'km'
    if css:
        cfbs = cfs+'_'+cdt0+'t0_'+cdte+'_'+css+'-'+crkm+'km'        
    print('    * which is original record '+str(krec)+' => date =',cdts,'=>',cfbs)
    cf_Q = './npz/Q-mesh_'+cfbs+'.npz'
    #
    return kr, krec, cdts, cdte, cfbs, cf_Q

#kr, krec, cdts, cdte, cfbs, cf_Q = _init4rec_( kr, pRec, pdate, cdt0, cfs, crkm, css )



def _QuadStat_( kr, QD ):

    # Some info about the spatial scales of quadrangles:
    print('\n *** About quadrangles at record #'+str(kr+1)+':')
    zsides = QD.lengths()
    zareas = QD.area()

    if np.any(zareas<0):
        print('ERROR: there are quads with a negative area :(')
        (idxFU,) = np.where( zareas<0 )
        print('  ==> '+str(len(idxFU))+' of such quads!')
        exit(0)
    
    if np.isnan(zareas).any():
        print('ERROR: there are "NANs" in `zarea`!!!!'); exit(0)
    zscale = np.sqrt(zareas)
    rl_average_side = np.mean(zsides)
    rl_average_scal = np.mean(zscale)
    rl_stdev_scal   = mjt.StdDev(rl_average_scal, zscale)
    rl_average_area = np.mean(zareas)
    print('  Quads @ rec '+str(kr)+' ==> scale: mean, StDev, min, max =',round(rl_average_scal,3), round(rl_stdev_scal,2),
          round(np.min(zscale),2), round(np.max(zscale),3),' km')
    print('  Quads @ rec '+str(kr)+' ==> average side length is '+str(round(rl_average_side,3))+' km')
    print('  Quads @ rec '+str(kr)+' ==> average area is '+str(round(rl_average_area,1))+' km^2')
    del zareas, zsides, zscale
    zdev = abs(rl_average_scal-reskm)
    if zdev > cfg.rc_tolQuadA:
        print(' ERROR: the mean scale is too different from the '+creskm
              +'km expected!!! (dev.=',round(zdev,2),' tol. = '+str(cfg.rc_tolQuadA)+'km)')
        exit(0)
    #
    return 0

    


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
    if Nrec!=2:
        print('ERROR: we expect only 2 records!!! Nrec =',Nrec)
        exit(0)    
    vRec = np.array( [ int(vcrec[i]) for i in range(Nrec) ], dtype=int )
    del vcrec

    mjt.chck4f(cf_nc_in)
    
    if not path.exists('./npz'): mkdir('./npz')
    if iplot>0:
        fdir = './figs/quadgener/'+str(reskm)+'km'
        makedirs( fdir, exist_ok=True )
    
    # Loading the data for the 2 selected records:
    Nt, NbP, corigin, lTimePos = mjt.GetDimNCdataMJT( cf_nc_in )

    if NbP<4:
        print('  \n *** Only '+str(NbP)+' points in the cloud! => exiting!!!');  exit(0)    
    
    if lTimePos: print(' *** "time_pos" data present in input netCDF file! => will use it!')

    kl = cfg.initialize( mode=quality_mode )
    kk = cfg.updateConfig4Scale( reskm, mode=quality_mode )

    zA_ideal = reskm*reskm ; # ideal Area of quadrangles to build!
    
    #print('\n *** Allowed deviation from '+creskm+' km for the MEAN scale of constructed quads (i.e. `sqrt(mean(Quad_areas))`) = ',cfg.rc_tolQuadA,'km')
    print('\n *** Max time deviation accepted across vertices of a Quad: `rc_t_dev_cancel` =',cfg.rc_t_dev_cancel,'s')
    print('\n *** Will only retain quadrangles with an area comprised between '
          +str(round(cfg.rc_Qarea_min,1))+' km^2 and '+str(round(cfg.rc_Qarea_max,1))+' km^2, Nominal='+str(zA_ideal)+'km^2')
    print('       + min and max angles, and ratio dev:', cfg.rc_Qang_min, cfg.rc_Qang_max, cfg.rc_dRatio_max)
    
    # We need a string to name the output files:
    kdt = -4
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
        print( ' * jrec = ',jrec, ', mean date =',e2c(vdate[jr]))    
        if jr==0:
            vIDs[:] = zIDs[:]
        else:
            if np.sum(zIDs[:]-vIDs[:])!=0:
                print('ERROR: ID fuck up in input file!') ; exit(0)
        jr=jr+1

        
    # Need some calendar info:
    NbDays = int( (vdate[1] - vdate[0]) / cfg.rc_day2sec )
    cdt1 = e2c(vdate[0] )
    cdt2 = e2c(vdate[-1])

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
        for jp in range(NbP):
            ztim[jp,:] = vdate[:]

    if idebug>0:
        for jr in range(Nrec):
            print('\n  DEBUG *** Record #'+str(vRec[jr])+':')
            for jc in range(NbP):
                print( ' * #'+str(jc)+' => Name: "'+zPnm[jc]+'": ID=',vIDs[jc],', X=',zXY[jc,0,jr],', Y=',zXY[jc,1,jr],
                       ', lon=',zGC[jc,0,jr],', lat=',zGC[jc,1,jr], ', time=',e2c(ztim[jc,jr]) )
                if str(vIDs[jc])!=zPnm[jc]:
                    print(' Fuck Up!!!! => vIDs[jc], zPnm[jc] =',vIDs[jc], zPnm[jc] ); exit(0)
        print('')

        
    cdate0 = _formDate_( vdate[0] )
    
    ##############
    # 1st record #
    ##############

    jr, jrec, cdats, cdate, cfbase1, cf_npzQ1 = _init4rec_( 0, vRec, vdate, cdate0, cfstr, creskm, crd_ss )
    
    #jr   = 0
    #jrec = vRec[jr]
    #print('\n\n *** QUAD-MESH GENERATION => record #'+str(jrec)+': record = '+str(jr))
    #cdats = e2c(vdate[jr])
    #cdate = _formDate_( vdate[jr] )
    #cfbase1 = cfstr+'_'+cdate0+'t0_'+cdate+'_'+creskm+'km'
    #if crd_ss:
    #    cfbase1 = cfstr+'_'+cdate0+'t0_'+cdate+'_'+crd_ss+'-'+creskm+'km'        
    #print('    * which is original record '+str(jrec)+' => date =',cdats,'=>',cfbase1)
    #cf_npzQ1 = './npz/Q-mesh_'+cfbase1+'.npz'
    
    print('\n *** Delaunay triangulation for 1st record!')

    cf_npzT1 = './npz/T-mesh_'+cfbase1+'.npz'

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
    xPxy1, vPids1, vTime1, xQpnts1, vQnam1 = mjt.Tri2Quad( TRIAS, anglRtri=(cfg.rc_Tang_min,cfg.rc_Tang_max),
                                                      ratioD=cfg.rc_dRatio_max, anglR=(cfg.rc_Qang_min,cfg.rc_Qang_max),
                                                      areaR=(cfg.rc_Qarea_min,cfg.rc_Qarea_max), areaIdeal=zA_ideal,
                                                      idbglev=idebug )

    lquit = ( np.shape(xPxy1)==(0,) or np.shape(xQpnts1)==(0,) )
    if not lquit:
        (nbP,_), (nbQ,_) = np.shape(xPxy1), np.shape(xQpnts1)                
        lquit = (nbQ<=0)
    if lquit:
        print('  \n *** 0 quads from triangles!!! => exiting!!!');  exit(0)
        
    print('\n *** After `Tri2Quad()`: we have '+str(nbQ)+' quadrangles out of '+str(nbP)+' points!\n')


    if corigin == 'RGPS':
        nbQo,nbPo = nbQ,nbP
        # Check if time deviation accros the 4 vertices of a quadrangle is too large, and cancel it if so...
        xPxy1, vPids1, vTime1, xQpnts1, vQnam1, _, _ = mjt.CancelQuadVrtcTDev( xPxy1, vPids1, vTime1, xQpnts1, vQnam1,
                                                                   cfg.rc_t_dev_cancel, mode=mode_ctl_vrtc )
        (nbP,_), (nbQ,_) = np.shape(xPxy1), np.shape(xQpnts1)
        print(' *** '+str(nbQo-nbQ)+' quads ('+str(nbPo-nbP)+' points) were supressed for inconsistent time across vertices! (t_dev<'
              +str(round(cfg.rc_t_dev_cancel/3600.,3))+'h)')
        
    # Save the triangular mesh info:
    mjt.SaveClassPolygon( cf_npzT1, TRIAS, ctype='T', origin=corigin )

    # To be used for other record, indices of Points to keep for Quads:
    _,ind2keep,_ = np.intersect1d(vIDs, vPids1, return_indices=True); # retain only indices of `vIDs` that exist in `vPids1`


    # Plots specific to `jrec==0`:
    if iplot>0:
        # We need to find a descent X and Y range for the figures:
        vrngX = mjt.roundAxisRange( TRIAS.PointXY[:,0], rndKM=50. )
        vrngY = mjt.roundAxisRange( TRIAS.PointXY[:,1], rndKM=50. )
        #
        if iplot>2:
            # Show triangles on a map:
            print('\n *** Launching Triangle plot!')
            kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig=fdir+'/fig01_Mesh_Triangles_'+cfbase1+'.png',
                                 ppntIDs=TRIAS.PointIDs,
                                 TriMesh=TRIAS.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig,
                                 rangeX=vrngX, rangeY=vrngY )
            

    #
    #if lExportCloudPoints:  
    #    # To be saved in a NC file:
    #    zpid  = np.zeros(nbP, dtype=int)
    #    zPxy  = np.zeros((nbP,2,Nrec))
    #    zPgc  = np.zeros((nbP,2,Nrec))
    #    zmask = np.zeros((nbP,Nrec))
    #    ztime = np.zeros((nbP,Nrec))
    #    zPxy[:,:,jr] = xPxy1[:,:]
    #    for jb in range(nbP):
    #        [zx,zy] = xPxy1[jb,:]
    #        #print(' zx,zy =',zx,zy)
    #        (idx,) = np.where( (np.round(zXY[:,0,jr],6) == round(zx,6)) & (np.round(zXY[:,1,jr],6) == round(zy,6)) )
    #        #print('idx =',idx)
    #        if len(idx)!=1:
    #            print('ERROR: `len(idx)!=1` !'); exit(0)
    #        zPgc[jb,:,jr] =  zGC[idx[0],:,jr]
    #        zpid[jb]      = vIDs[idx[0]]
    #        zmask[jb,jr] = zmsk[idx[0],jr]
    #        ztime[jb,jr] =  ztim[idx[0],jr]
        
    
    
    #######################################################################################################

    ##############
    # 2nd record #
    ##############

    jr, jrec, cdats, cdate, cfbase2, cf_npzQ2 = _init4rec_( 1, vRec, vdate, cdate0, cfstr, creskm, crd_ss )

    print('\n *** Quad recycling for jr=',jr)

    QUADS1 = mjt.Quadrangle( xPxy1, xQpnts1, vPids1, vTime1, vQnam1, date=cdats, origin=corigin, reskm_nmnl=reskm )

    if iplot>1:
        print('\n *** Launching Quad-only plot for Quads at 1st record (BEFORE 2nd pass REMOVAL )!')
        kk = mjt.ShowTQMesh( QUADS1.PointXY[:,0], QUADS1.PointXY[:,1], cfig=fdir+'/fig03_BEFORE-clean_QUADS1_'+cfbase1+'.png',
                             ppntIDs=QUADS1.PointIDs, QuadMesh=QUADS1.MeshVrtcPntIdx, qIDs=QUADS1.QuadIDs,
                             lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
    
    # Recycling Quads found at 1st record (QUADS1):
    xPxy2, vTime2, xQpnts2, vPids2, vQnam2, vQIDs2  = mjt.RecycleQuads( zXY[:,:,jr], ztim[:,jr], vIDs, QUADS1,  iverbose=idebug )

    # Time to spot possible issues in these recycled quads at 2nd record and fix them...
    ifix = 0 ; # if ifix>0 later, then QUADS1 has to be rebuilt!
    
    # Fix #1: possible that some quads that were ok at the first record have become twisted (with a negative area) now
    # at the second record:
    nbQo,nbPo = nbQ,nbP
    zQcoor = np.array( [ xPxy2[xQpnts2[jQ,:],:] for jQ in range(nbQ) ])
    zareas = mjt.QuadsAreas( zQcoor )
    if np.any(zareas < 0.):
        print('WARNING: some quads have become twisted (negative area) in the second record (recycling)!')
        (idxFU,) = np.where(zareas < 0.)
        nFU = len(idxFU)
        print('       => '+str(nFU)+' of such quads...')
        print('       => calling `CancelQuadNegArea()` !')
        #
        xPxy2, vPids2, vTime2, xQpnts2, vQnam2, idxQk, idxPk = mjt.CancelQuadNegArea( xPxy2, vPids2, vTime2, xQpnts2, vQnam2 )
        #
        (nbP,_), (nbQ,_) = np.shape(xPxy2), np.shape(xQpnts2)
        if nbQo-nbQ != nFU:
            print('ERROR: `nbQo-nbQ != nFU` (Areas) !'); exit(0)
        #
        print('    * '+str(nbQo-nbQ)+' quads in second records had to be supressed because of negative area!!!')
        print('       => cancelling the same quads at 1st record as well!')
        xPxy1, vPids1, vTime1, xQpnts1, vQnam1, _ = mjt.KeepSpcfdQuads( idxQk, xPxy1, vPids1, vTime1, xQpnts1, vQnam1 )
        ifix+=1
    #
    del zQcoor, zareas
        
    # Fix #2: here too, we must ensure that no quadrangle bears too much time inconsitency across vertices:
    if corigin == 'RGPS':
        nbQo,nbPo = nbQ,nbP
        # Check if time deviation accros the 4 vertices of a quadrangle is too large, and cancel it if so...
        xPxy2, vPids2, vTime2, xQpnts2, vQnam2, idxQk, idxPk = mjt.CancelQuadVrtcTDev( xPxy2, vPids2, vTime2, xQpnts2, vQnam2,
                                                                   cfg.rc_t_dev_cancel, mode=mode_ctl_vrtc )
        (nbP,_), (nbQ,_) = np.shape(xPxy2), np.shape(xQpnts2)
        if nbQ<nbQo:
            print('    * '+str(nbQo-nbQ)+' quads in second records had to be supressed because of vertices time inconsistencies!!!')
            print('       => cancelling the same quads at 1st record as well!')
            xPxy1, vPids1, vTime1, xQpnts1, vQnam1, _ = mjt.KeepSpcfdQuads( idxQk, xPxy1, vPids1, vTime1, xQpnts1, vQnam1 )
            ifix+=1


    
    # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
    if ifix>0:
        QUADS1 = mjt.Quadrangle( xPxy1, xQpnts1, vPids1, vTime1, vQnam1, date=cdats, origin=corigin, reskm_nmnl=reskm )        
    QUADS2     = mjt.Quadrangle( xPxy2, xQpnts2, vPids2, vTime2, vQnam2, vQIDs=vQIDs2, date=cdats, origin=corigin, reskm_nmnl=reskm )

    del xPxy1, xQpnts1, vPids1, vTime1, vQnam1
    del xPxy2, xQpnts2, vPids2, vTime2, vQnam2, vQIDs2

    
    k1 = _QuadStat_( 0, QUADS1 )
    k2 = _QuadStat_( 1, QUADS2 )
    
    # Save the quadrangular mesh info:
    mjt.SaveClassPolygon( cf_npzQ1, QUADS1, ctype='Q', origin=corigin, reskm_nmnl=reskm )
    mjt.SaveClassPolygon( cf_npzQ2, QUADS2, ctype='Q', origin=corigin, reskm_nmnl=reskm )
    
        
    if iplot>0:

        if iplot>1:
            # Show triangles together with the quadrangles on a map:
            print('\n *** Launching Triangle+Quad plot!')
            kk = mjt.ShowTQMesh( TRIAS.PointXY[:,0], TRIAS.PointXY[:,1], cfig=fdir+'/fig02_TRI1_'+cfbase1+'.png',
                                 ppntIDs=TRIAS.PointIDs, TriMesh=TRIAS.MeshVrtcPntIdx,
                                 pX_Q=QUADS1.PointXY[:,0], pY_Q=QUADS1.PointXY[:,1], QuadMesh=QUADS1.MeshVrtcPntIdx,
                                 lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )
            
        if iplot>2:
            # Show only points composing the quadrangles:
            kk = mjt.ShowTQMesh( QUADS1.PointXY[:,0], QUADS1.PointXY[:,1], cfig=fdir+'/fig03a_QUADS1_'+cfbase1+'.png',
                                 ppntIDs=QUADS1.PointIDs, lGeoCoor=False, zoom=rzoom_fig )
        
        # Show only the quads with only the points that define them:
        print('\n *** Launching Quad-only plot for Quads at 1st record!')
        kk = mjt.ShowTQMesh( QUADS1.PointXY[:,0], QUADS1.PointXY[:,1], cfig=fdir+'/fig03_QUADS1_'+cfbase1+'.png',
                             ppntIDs=QUADS1.PointIDs, QuadMesh=QUADS1.MeshVrtcPntIdx, qIDs=QUADS1.QuadIDs,
                             lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )

        # 2nd record:
        print('\n *** Launching Quad-only plot for Quads at 2nd record!')
        kk = mjt.ShowTQMesh( QUADS2.PointXY[:,0], QUADS2.PointXY[:,1], cfig=fdir+'/fig03_QUADS2_'+cfbase2+'.png',
                             ppntIDs=QUADS2.PointIDs, QuadMesh=QUADS2.MeshVrtcPntIdx, qIDs=QUADS2.QuadIDs,
                             lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )

        
        #kk = mjt.ShowQuads( QUA.PointXY[QUA.MeshVrtcPntIdx,:], cfig=fdir+'/fig04_Quadrangles_'+cfbase+'.png', rangeX=vrngX, rangeY=vrngY )



    del TRIAS, QUADS1, QUADS2

    print('\n --- Over!\n')
    
    exit(0)
    

    #if lExportCloudPoints:
    #    # To be save in a NC file:
    #    zPxy[:,:,jr] = xPxy2[:,:]
    #    for jb in range(nbP):
    #        [zx,zy] = xPxy2[jb,:]
    #        (idx,) = np.where( (np.round(zXY[:,0,jr],6) == round(zx,6)) & (np.round(zXY[:,1,jr],6) == round(zy,6)) )
    #        if len(idx)!=1:
    #            print('ERROR: `len(idx)!=1` !'); exit(0)
    #        zPgc[jb,:,jr] =  zGC[idx[0],:,jr]
    #        zpid[jb]      = vIDs[idx[0]]
    #        zmask[jb,jr]  = zmsk[idx[0],jr]
    #        ztime[jb,jr]  =  ztim[idx[0],jr]


    #if lExportCloudPoints:
    #    cf_nc_out = str.replace( cf_nc_in, '.nc', '_postQG.nc' )
    #
    #    zYkm = np.zeros( (Nrec,nbP) )
    #    zXkm = np.zeros( (Nrec,nbP) )
    #    zlat = np.zeros( (Nrec,nbP) )
    #    zlon = np.zeros( (Nrec,nbP) )
    #
    #    zYkm[:,:] = zPxy[:,1,:].T
    #    zXkm[:,:] = zPxy[:,0,:].T
    #    zlat[:,:] = zPgc[:,1,:].T
    #    zlon[:,:] = zPgc[:,0,:].T
    #
    #    print('   * saving '+cf_nc_out)
    #
    #    kk = mjt.ncSaveCloudBuoys( cf_nc_out, vdate, zpid, zYkm, zXkm, zlat, zlon, mask=zmask.T,
    #                               xtime=ztime.T, fillVal=mjt.FillValue, corigin=corigin )
    #    if iplot>0:
    #        # Ploting what we saved in NC file:
    #        for jr in range(Nrec):
    #            cfig = str.replace( cf_nc_out, '.nc', '_postQG'+'_rec%3.3i'%(jr)+'.png' )
    #            cfig = str.replace( cfig , './nc', fdir )
    #            kk = mjt.ShowBuoysMap( vdate[jr], zlon[jr,:], zlat[jr,:], pvIDs=vIDs[:], cfig=cfig,
    #                                   cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )
    #

    

    

                
        
    
    




    
