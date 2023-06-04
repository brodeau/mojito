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
    cdate = str.replace( e2c(idate, precision='D'), '-', '')
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




if __name__ == '__main__':    
    
    crd_ss = None
    
    if not len(argv) in [6,7]:
        print('Usage: '+argv[0]+' <file_mojito.nc> <records to use (C), comma-separated> <npz_model_Quad_file> <basic_res_km> <mode (rgps,model,xlose)> (<rd_ss_km>)')
        exit(0)

    cf_nc_in = argv[1]
    lstrec  = argv[2]
    cf_Q_npz = argv[3]
    creskm = argv[4]
    reskm = int(creskm)
    quality_mode = argv[5]
    ik = cfg.controlModeName( path.basename(__file__), quality_mode )

    if len(argv)==7:
        crd_ss = argv[6]

    
    vcrec = split(',',lstrec)
    Nrec  = len(vcrec)
    if Nrec!=2:
        print('ERROR: we expect only 2 records!!! Nrec =',Nrec)
        exit(0)    
    vRec = np.array( [ int(vcrec[i]) for i in range(Nrec) ], dtype=int )
    del vcrec

    mjt.chck4f(cf_nc_in)
    mjt.chck4f(cf_Q_npz)
    
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

    print('\n *** Quad recycling for 1st record jr=',jr)
    print('      => based on Quads specs found in file: '+cf_Q_npz)

    QUADStmplt = mjt.LoadClassPolygon( cf_Q_npz, ctype='Q' )


    nbP = QUADStmplt.nP
    nbQ = QUADStmplt.nQ
    
    print('      => number of Quads and Points in this template npz file:',nbQ,nbP)
    
    # Recycling Quads of template file using new info read in netCDF file:
    xPxy1, vTime1, xQpnts1, vPids1, vQnam1, vQIDs1  = mjt.RecycleQuads( zXY[:,:,jr], ztim[:,jr], vIDs, QUADStmplt,  iverbose=idebug )
    print('      => Recycling done for 1st record! :)\n')
    del QUADStmplt

    QUADS1 = mjt.Quadrangle( xPxy1, xQpnts1, vPids1, vTime1, vQnam1, date=cdats, origin=corigin, reskm_nmnl=reskm )        
    print('      => "QUADS1" constructed! :)\n')
    
    if iplot>0:
        vrngX = mjt.roundAxisRange( QUADS1.PointXY[:,0], rndKM=50. )
        vrngY = mjt.roundAxisRange( QUADS1.PointXY[:,1], rndKM=50. )
    if iplot>1:
        print('\n *** Launching Quad-only plot for Quads at 1st record (BEFORE 2nd pass REMOVAL )!')
        kk = mjt.ShowTQMesh( QUADS1.PointXY[:,0], QUADS1.PointXY[:,1], cfig=fdir+'/fig03_BEFORE-clean_QUADS1_'+cfbase1+'.png',
                             ppntIDs=QUADS1.PointIDs, QuadMesh=QUADS1.MeshVrtcPntIdx, qIDs=QUADS1.QuadIDs,
                             lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY )

    k1 = mjt.QuadStat( 0, QUADS1, resolkm=reskm, tolArea=cfg.rc_tolQuadA )
        
    #######################################################################################################

    ##############
    # 2nd record #
    ##############    
    jr, jrec, cdats, cdate, cfbase2, cf_npzQ2 = _init4rec_( 1, vRec, vdate, cdate0, cfstr, creskm, crd_ss )

    print('\n *** Quad recycling for jr=',jr)

    
    # Recycling Quads found at 1st record (QUADS1):
    xPxy2, vTime2, xQpnts2, vPids2, vQnam2, vQIDs2  = mjt.RecycleQuads( zXY[:,:,jr], ztim[:,jr], vIDs, QUADS1,  iverbose=idebug )
    print('      => Recycling done for 2nd record! :)\n')
    
    # Time to spot possible issues in these recycled quads at 2nd record and fix them...
    ifix = 0 ; # if ifix>0 later, then QUADS1 has to be rebuilt!
    
    # Fix #1: possible that some quads that were ok at the first record have become twisted (with a negative area) now
    # at the second record:
    for ii in range(2):
        nbQo,nbPo = nbQ,nbP
        zQcoor = np.array( [ xPxy2[xQpnts2[jQ,:],:] for jQ in range(nbQ) ])
        zareas = mjt.QuadsAreas( zQcoor )
        if np.any(zareas <= 0.):
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
        
    
    # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
    if ifix>0:
        print('WARNING: fixing (updating) QUADS1 !!!')
        QUADS1 = mjt.Quadrangle( xPxy1, xQpnts1, vPids1, vTime1, vQnam1, date=cdats, origin=corigin, reskm_nmnl=reskm )
        k1 = mjt.QuadStat( 0, QUADS1, resolkm=reskm, tolArea=cfg.rc_tolQuadA )
        
    QUADS2     = mjt.Quadrangle( xPxy2, xQpnts2, vPids2, vTime2, vQnam2, vQIDs=vQIDs2, date=cdats, origin=corigin, reskm_nmnl=reskm )
    print('      => "QUADS2" constructed! :)\n')
    k2 = mjt.QuadStat( 1, QUADS2, resolkm=reskm, tolArea=cfg.rc_tolQuadA )
    
    del xPxy1, xQpnts1, vPids1, vTime1, vQnam1
    del xPxy2, xQpnts2, vPids2, vTime2, vQnam2, vQIDs2
    
    # Save the quadrangular mesh info:
    mjt.SaveClassPolygon( cf_npzQ1, QUADS1, ctype='Q', origin=corigin, reskm_nmnl=reskm )
    mjt.SaveClassPolygon( cf_npzQ2, QUADS2, ctype='Q', origin=corigin, reskm_nmnl=reskm )
    
        
    if iplot>0:

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



    del QUADS1, QUADS2

    print('\n --- Over!\n')
    

    
    




    
