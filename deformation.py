#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
    # The present script should only compute deformations, and nothing else!
    # As such, it is not expected to do any selection of what deformations 
    # are acceptable or not! => should be done elsewhere

'''

from sys import argv, exit
from os import path, mkdir, makedirs
import numpy as np
from re import split

from scipy.spatial import Delaunay

import mojito   as mjt
from mojito import config as cfg
from mojito.util import epoch2clock as e2c

idebug=0
iplot=1 ; NameArcticProj='SmallArctic'
zoom=2

lExportNC4Seed = False ; # Triggers the saving of the "mojito" netCDF and npz files that can be used to seed in `sitrack` the
#                       # points that define the quads for which a deformation was computed, and save a npz file containaing the
#                       # Quad class to allow for reconstruction of the same quads...




if __name__ == '__main__':


    if not len(argv) in [5,6]:
        print('Usage: '+argv[0]+' <file_Q_mesh_N1.npz> <file_Q_mesh_N2.npz> <time_dev_from_mean_allowed (s)> <mode (rgps,model,xlose)> (<E:[export{ Quad info]>)')
        exit(0)
    cf_Q1 = argv[1]
    cf_Q2 = argv[2]
    time_dev_max= float(argv[3])
    quality_mode = argv[4]

    if len(argv)==6:
        lExportNC4Seed = (argv[5] in ['E','e'] )

    ik = cfg.controlModeName( path.basename(__file__), quality_mode )
    print('\n *** Max time_dev_from_mean_allowed =',time_dev_max/3600,'hours')

    # Reading the quad meshes in both npz files:
    QUA1 = mjt.LoadClassPolygon( cf_Q1, ctype='Q' )
    nP1,nQ1 = QUA1.nP,QUA1.nQ
    print('        =>  '+str(nQ1)+' Quads constructed on '+str(nP1)+' points.')
    
    QUA2 = mjt.LoadClassPolygon( cf_Q2, ctype='Q' )    
    nP2,nQ2 = QUA2.nP,QUA2.nQ
    print('        =>  '+str(nQ2)+' Quads constructed on '+str(nP2)+' points.')

    
    # Debug have a look at the times of all points and get the actual mean time for each file:
    rTm1, rStD1 = mjt.CheckTimeConsistencyQuads(1, QUA1, time_dev_max, iverbose=idebug)
    rTm2, rStD2 = mjt.CheckTimeConsistencyQuads(2, QUA2, time_dev_max, iverbose=idebug)

    itimeC = int( 0.5*(rTm1+rTm2) )
    ctimeC = e2c(itimeC)
    print('\n *** Deformations will be calculated at about: '+ctimeC+'\n')
    rdt = rTm2 - rTm1

    vclck = split('_',ctimeC)
    chh, cmm = split(':',vclck[1])[0], split(':',vclck[1])[1]
    cclck = e2c(itimeC, precision='D', frmt='nodash') ; #str.replace( vclck[0],'-','')+'-'+chh+'h'+cmm

    reskm, reskm2 = QUA1.reskm_nmnl, QUA2.reskm_nmnl
    if reskm != reskm2:
        mjt.printEE('quads do not have the same nominal resolution in the 2 files:',reskm, reskm2)        
    print('\n *** Nominal resolution for the Quads in both files =',reskm,'km')


    makedirs( './npz', exist_ok=True )    
    if iplot>0:
        cfdir = './figs/deformation/'+str(reskm)+'km'
        makedirs( cfdir, exist_ok=True )

    # Comprehensive name for npz and figs to save later on:
    corigin = QUA1.origin
    
    k1 = cfg.initialize(                mode=quality_mode )
    k2 = cfg.updateConfig4Scale( reskm, mode=quality_mode )
    print(' *** Min and max deformation allowed:',cfg.rc_tot_min, cfg.rc_tot_max,' days^-1 !')
    if not cfg.lc_accurate_time:
        print(' *** Time step to be used: `dt` = '+str(round(rdt,2))+' = '+str(round(rdt/cfg.rc_day2sec,2))+' days')

        
    # Some info from npz file name: #fixme: too much dependency on file name...
    cf1, cf2 = path.basename(cf_Q1), path.basename(cf_Q2)
    if corigin == 'RGPS':
        cbatch =  split('_',cf1)[2]
        cdtbin =  '_'+split('_',cf1)[3]
    elif split('_',corigin)[0] == 'NEMO-SI3':
        cbatch =  split('_',cf1)[-5]
        cdtbin = '_'+split('_',cf1)[-4]
    if not cdtbin[1:3] in ['dt','No']:
        mjt.printEE('we could not figure out `cdtbin`!')


    # Test for coarsening realisations:
    cr1, cr2 = '', ''
    cc1, cc2 = split('-',split('_',cf1)[-1]), split('-',split('_',cf2)[-1])
    if len(cc1)==2: cr1 = cc1[0]+'-'
    if len(cc2)==2: cr2 = cc2[0]+'-'
    if cr1!='' and cr1!=cr2:
        mjt.printEE('with the realisation coarsening res:',cr1,cr2)
    else:
        print('     => realisation with `rd_ss` =',cc1[0],'km !!!')

    cfnm0 = corigin+'_'+cbatch+cdtbin
    cfnm  = cfnm0+'_'+cclck
    cfnm += '_'+cr1+str(reskm)+'km'


    if  cdtbin=='_NoBin':
        dtbin = 0
    else:
        dtbin = int(cdtbin[3:])*3600
        print('\n *** width of time bin used in RGPS =',dtbin/3600,'hours!')

    if cfg.lc_accurate_time:
        figSfx='_tbuoy.png'
    else:
        figSfx='_tglob.png'

    
    vnm, idxK1, idxK2 = np.intersect1d( QUA1.QuadNames, QUA2.QuadNames, assume_unique=True, return_indices=True )
        
    if nQ1 == nQ2:
        print('\n *** Great! Both files have the same number of Quads!!!')
        nQ = nQ1
    else:
        print('\n [WARNING]: the 2 files/records do not have the same number of Quads!!!')
        # The Quads we retain, i.e. those who exist in both snapshots:        
        nQ = len(vnm) ; # also = len(vidx*)
        znm, zidx1, zidx2 = np.intersect1d( QUA1.QuadIDs, QUA2.QuadIDs, assume_unique=True, return_indices=True )
        nQ2 = len(znm)
        if nQ!=nQ2 or np.sum(zidx1-idxK1)!=0 or np.sum(zidx2-idxK2)!=0:
            mjt.printEE('we do not get the same info based on Quad names and Quad Ids !!!')
        #
    print('       => there are '+str(nQ)+' Quads common to the 2 records!\n')


    # Now for some weird reasons time of a given buoy can be the same in the 2 quads:
    zTime1 = QUA1.MeshVrtcPntTime()[idxK1,:]
    zTime2 = QUA2.MeshVrtcPntTime()[idxK2,:]
    zdT = zTime2-zTime1
    if np.any(zdT==0.):
        mjt.printEE('time for some buoys is the same in the 2 records!')
        # => this should be fixed at the quad generation level!!! No here!!!
    del zTime1, zTime2, zdT


    # Okay, no we can start the real shit...
    
    # Coordinates of the 4 points of quadrangles for the 2 consecutive records:
    zXY1 = QUA1.MeshPointXY[idxK1,:,:].copy() ; #  km !
    zXY2 = QUA2.MeshPointXY[idxK2,:,:].copy() ; #  km !

    # Computation of partial derivative of velocity vector constructed from the 2 consecutive positions:
    if cfg.lc_accurate_time:
        # Time of poins of the 4 points of quadrangles for the 2 consecutive records:
        zTime1 = QUA1.MeshVrtcPntTime()[idxK1,:]
        zTime2 = QUA2.MeshVrtcPntTime()[idxK2,:]
        #
        zX, zY, zU, zV, zdUdxy, zdVdxy, zAq = mjt.PDVfromPos( 1, zXY1, zXY2, QUA1.area()[idxK1], QUA2.area()[idxK2],
                                                              xtime1=zTime1, xtime2=zTime2, iverbose=idebug )
        # => having rdt=1 will yield fuck-up fields if used, and must not be used!

    else:

        zX, zY, zU, zV, zdUdxy, zdVdxy, zAq = mjt.PDVfromPos( rdt, zXY1, zXY2, QUA1.area()[idxK1], QUA2.area()[idxK2],
                                                              iverbose=idebug )

    del zXY1, zXY2
    
    # zX, zY => positions  of the 4 vertices at center of time interval!
    # zU, zV => velocities of the 4 vertices at center of time interval!
    #   zAq  => Quad area used in the estimation of `zdUdxy, zdVdxy` !

    # For plots:
    zrx = [ np.min(zX)-25. , np.max(zX)+25. ]
    zry = [ np.min(zY)-25. , np.max(zY)+25. ]


    # Coordinates of barycenter of Quads at center of time interval:
    zXc = np.mean( zX[:,:], axis=1 )
    zYc = np.mean( zY[:,:], axis=1 )


    if idebug>1 and iplot>0:
        # Velocities of barycenter of Quads at center of time interval:
        zUc = np.mean( zU[:,:], axis=1 )
        zVc = np.mean( zV[:,:], axis=1 )
        #
        mjt.ShowDefQuad( zXc, zYc, zUc, cfig=cfdir+'/zvUc_'+cfnm+figSfx, cwhat='Uc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, unit='m/s' )
        mjt.ShowDefQuad( zXc, zYc, zVc, cfig=cfdir+'/zvVc_'+cfnm+figSfx, cwhat='Vc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, unit='m/s' )
        mjt.ShowDefQuad( zXc, zYc, np.sqrt(zUc*zUc+zVc*zVc), cfig=cfdir+'/zUMc_'+cfnm+figSfx, cwhat='UMc',
                             pFmin=0., pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, unit='m/s' )
        del zUc, zVc

    # Divergence aka SigmaI:
    zdiv = mjt.DivPDV(zdUdxy, zdVdxy)

    # `Maximum Shear Strain Rate`**2:
    zshr2 = mjt.Shear2PDV(zdUdxy, zdVdxy )

    del zdUdxy, zdVdxy

    (nD,) = np.shape(zdiv)
    if  nQ != nD:
        mjt.printEE('`Q != nD` !!!')
    if  np.shape(zshr2) != (nD,):
        mjt.printEE('`np.shape(zdiv) != np.shape(zshr2)` !!!')

    # Maximum Shear Strain Rate aka SigmaII:
    zshr = np.zeros(nD)
    zshr = np.sqrt(zshr2)

    # Total deformation rate:
    ztot = np.zeros(nD)
    ztot[:] = np.sqrt( zdiv[:]*zdiv[:] + zshr2[:] )
    del zshr2


    #idxKeep = np.arange(nD,dtype=int) ; # default we keep everything
    #print(idxKeep[::10])
    #exit(0);#lili

    lNeedClean = False

    
    # Non-realistic / error extreme values in computed deformation:
    if quality_mode in ['rgps','rgps_map','rgps_track']:
        ztotdm1 = ztot*cfg.rc_day2sec ; # same but in days^-1 !

        lNeedClean = np.any( (ztotdm1 < cfg.rc_tot_min) | (ztotdm1 > cfg.rc_tot_max) )
        
        if lNeedClean:
            # Must get rid of extremely small deformation (if RGPS! if not => `rc_div_min, rc_shr_min, rc_tot_min` taken ridiculously tiny!)
            (idxKeep,) = np.where( (ztotdm1 >= cfg.rc_tot_min) & (ztotdm1 <= cfg.rc_tot_max) )
            #
            nDn = len(idxKeep)
            if nDn < nD:
                print('\n *** EXCLUDING '+str(nD-nDn)+' deformation points (quads) because of excessively small deformation rate! (because mode='+quality_mode+')')
                zdiv = zdiv[idxKeep]
                zshr = zshr[idxKeep]
                ztot = ztot[idxKeep]
                zXc  =  zXc[idxKeep]
                zYc  =  zYc[idxKeep]
                zAq  =  zAq[idxKeep]
                nD   = nDn
                zX   = zX[idxKeep,:]
                zY   = zY[idxKeep,:]
        else:
            print('\n *** Great! No extreme deformation values were found!')
            idxKeep = np.arange(nD, dtype=int)
            
        del ztotdm1
        #
        print('        =>  '+str(nD)+' Quads left.')
    
    # Save the deformation data:
    np.savez_compressed( './npz/DEFORMATIONS_'+cfnm+'.npz', time=itimeC, date=ctimeC, Npoints=nD,
                         Xc=zXc, Yc=zYc, X4=zX, Y4=zY, divergence=zdiv, shear=zshr, total=ztot,
                         quadArea=zAq, origin=corigin, reskm_nmnl=reskm )


    # Save the "mojito" netCDF and npz Quad class files that can be used to seed in `sitrack` for only the quads that had a reasonable deformation
    # and reconstruct the same quads
    if lExportNC4Seed:

        # Last-man standing quad indices:
        idxQ1, idxQ2 = idxK1[idxKeep], idxK2[idxKeep]
        if np.shape(idxQ1)!=(nD,) or np.shape(idxQ2)!=(nD,):
            mjt.printEE('fuck-up #1')

        zPXY1, zPids1, zTime1, zQpnts1, zQnam1, _ = mjt.KeepSpcfdQuads( idxQ1, QUA1.PointXY, QUA1.PointIDs, QUA1.PointTime, QUA1.MeshVrtcPntIdx, QUA1.QuadNames )
        zPXY2, zPids2, zTime2, zQpnts2, zQnam2, _ = mjt.KeepSpcfdQuads( idxQ2, QUA2.PointXY, QUA2.PointIDs, QUA2.PointTime, QUA2.MeshVrtcPntIdx, QUA2.QuadNames )
        if any(zPids2-zPids1!=0):
            mjt.printEE('fuck-up #2')
        
        QR1 = mjt.Quadrangle( zPXY1, zQpnts1, zPids1, zTime1, zQnam1, origin='RGPS', reskm_nmnl=reskm )
        if idebug>0 and iplot>0:
            QR2 = mjt.Quadrangle( zPXY2, zQpnts2, zPids2, zTime2, zQnam2, origin='RGPS', reskm_nmnl=reskm )
                
        Nrec, Np, Nq = 2, len(zPids1), nD
        vtim = np.array([ np.mean(zTime1), np.mean(zTime2) ], dtype=int) ; # Fill `vtim` with mean date at each record:
        xtim = np.array([ zTime1, zTime2 ], dtype=int)        
        zPGC1, zPGC2 = mjt.CartNPSkm2Geo1D( zPXY1, convArray='F' ), mjt.CartNPSkm2Geo1D( zPXY2, convArray='F' )        
        zPXY, zPGC = np.array([ zPXY1, zPXY2]), np.array([ zPGC1, zPGC2])
        del zPGC1, zPGC2

        cdt1 = e2c(vtim[0], precision='D', frmt='nodash')
        cfnm1  = cfnm0+'_'+cdt1+'_'+cr1+str(reskm)+'km'
        cfnm2  = cfnm0+'_'+cdt1+'_'+e2c(vtim[1], precision='D', frmt='nodash')+'_'+cr1+str(reskm)+'km'

        k1 = mjt.QuadStat( 0, QR1, resolkm=reskm, tolArea=cfg.rc_maxDevMeanAreaQuads )        
        mjt.SaveClassPolygon( './npz/QUADSofDEF_'+cfnm1+'.npz', QR1, ctype='Q', origin='RGPS', reskm_nmnl=reskm )

        makedirs( './nc', exist_ok=True )    
        kk = mjt.ncSaveCloudBuoys( './nc/PointsOfQuadsOfDEF_'+cfnm2+'.nc', vtim, zPids1, zPXY[:,:,1], zPXY[:,:,0], zPGC[:,:,1], zPGC[:,:,0],
                                   xtime=xtim, fillVal=mjt.FillValue, corigin='RGPS' )
        del zPGC, xtim
        print('\n *** All info necessary to seed the points and reconstruct the Quads involved in computed deformation cells is saved!')
        print('         => that is '+str(Nq)+' Quads constructed on '+str(Np)+' points.\n')

    
    # Some plots:
    if iplot>0:

        if corigin != 'RGPS':
            corigin = str.replace( corigin,'NEMO-','')
            corigin = str.replace( corigin,'_NANUK4_','-')

        cresinfo = '('+str(reskm)+' km)'

        nmproj=NameArcticProj

        # Filled quads projected on the Arctic map:
        mjt.ShowDefQuadGeoArctic( zX, zY, cfg.rc_day2sec*zdiv, cfig=cfdir+'/map_zd_'+cfnm+'_Divergence'+figSfx, nmproj=NameArcticProj, cwhat='div',
                                  pFmin=-cfg.rc_div_max_fig, pFmax=cfg.rc_div_max_fig, pdF=cfg.rc_df_fig, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': divergence '+cresinfo, idate=itimeC, edgecolor=None )

        mjt.ShowDefQuadGeoArctic( zX, zY, cfg.rc_day2sec*zshr, cfig=cfdir+'/map_zs_'+cfnm+'_Shear'+figSfx,      nmproj=NameArcticProj, cwhat='shr',
                                  pFmin=0.,      pFmax=cfg.rc_shr_max_fig, pdF=cfg.rc_df_fig, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': shear '+cresinfo, idate=itimeC, edgecolor=None )


        cix = ''
        cix ='max%5.5i'%(int( round( 10000.*np.max(cfg.rc_day2sec*ztot) , 0) ))+'_'

        mjt.ShowDefQuadGeoArctic( zX, zY, cfg.rc_day2sec*ztot, cfig=cfdir+'/map_zt_'+cix+cfnm+'_Total'+figSfx,      nmproj=NameArcticProj, cwhat='tot',
                                  pFmin=0.,      pFmax=cfg.rc_tot_max_fig, pdF=cfg.rc_df_fig, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': total deformation '+cresinfo, idate=itimeC, edgecolor=None )


        if lExportNC4Seed and idebug>0:
            vrngX = mjt.roundAxisRange( QR1.PointXY[:,0], rndKM=50. )
            vrngY = mjt.roundAxisRange( QR1.PointXY[:,1], rndKM=50. )
            #
            jr=0
            kf = mjt.ShowTQMesh( zPXY[jr,:,0], zPXY[jr,:,1], cfig=cfdir+'/RECONSTRUCTED_remaining_Quads_'+cfnm2+'_rec%2.2i'%(jr)+'.png',
                                 ppntIDs=zPids1, QuadMesh=QR1.MeshVrtcPntIdx, qIDs=QR1.QuadIDs,
                                 lGeoCoor=False, zoom=4, rangeX=vrngX, rangeY=vrngY, lShowIDs=True )
            jr=1
            kf = mjt.ShowTQMesh( zPXY[jr,:,0], zPXY[jr,:,1], cfig=cfdir+'/RECONSTRUCTED_remaining_Quads_'+cfnm2+'_rec%2.2i'%(jr)+'.png',
                                 ppntIDs=zPids2, QuadMesh=QR2.MeshVrtcPntIdx, qIDs=QR2.QuadIDs,
                                 lGeoCoor=False, zoom=4, rangeX=vrngX, rangeY=vrngY, lShowIDs=True )
            


        
    if iplot>1:
        # Filled quads projected on RGPS projection (Cartesian):
        mjt.ShowDefQuad( zX, zY, cfg.rc_day2sec*zdiv, cfig=cfdir+'/zd_'+cfnm+'_Divergence'+figSfx, cwhat='div',
                                  pFmin=-cfg.rc_div_max_fig, pFmax=cfg.rc_div_max_fig, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': divergence '+cresinfo )

        mjt.ShowDefQuad( zX, zY, cfg.rc_day2sec*zshr, cfig=cfdir+'/zs_'+cfnm+'_Shear'+figSfx,      cwhat='shr',
                                  pFmin=0.,      pFmax=cfg.rc_shr_max_fig,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': shear '+cresinfo )

        mjt.ShowDefQuad( zX, zY, cfg.rc_day2sec*zshr, cfig=cfdir+'/zt_'+cfnm+'_Total'+figSfx,      cwhat='tot',
                                  pFmin=0.,      pFmax=cfg.rc_tot_max_fig,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': total deformation '+cresinfo )


        ###
