#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

from scipy.spatial import Delaunay

import mojito   as mjt
from mojito import config as cfg

idebug=0

iplot=1 ; NameArcticProj='SmallArctic'
zoom=2

if __name__ == '__main__':

    kk = cfg.initialize()
    
    if not len(argv) in [4]:
        print('Usage: '+argv[0]+' <file_Q_mesh_N1.npz> <file_Q_mesh_N2.npz> <time_dev_from_mean_allowed (s)>')
        exit(0)
    cf_Q1 = argv[1]
    cf_Q2 = argv[2]
    time_dev_max= float(argv[3])

    if iplot>0:
        cdir = './figs/deformation'
        for cc in ['./figs',cdir]:
            if not path.exists(cc): mkdir(cc)
        
    print('\n *** Max time_dev_from_mean_allowed =',time_dev_max/3600,'hours')
        
    # Reading the quad meshes in both npz files:
    QUA1 = mjt.LoadClassPolygon( cf_Q1, ctype='Q' )
    QUA2 = mjt.LoadClassPolygon( cf_Q2, ctype='Q' )

    # Debug have a look at the times of all points and get the actual mean time for each file:
    rTm1, rStD1 = mjt.CheckTimeConsistencyQuads(1, QUA1, time_dev_max, iverbose=idebug)
    rTm2, rStD2 = mjt.CheckTimeConsistencyQuads(2, QUA2, time_dev_max, iverbose=idebug)
    
    rtimeC = 0.5*(rTm1+rTm2)
    ctimeC = mjt.epoch2clock(rtimeC)
    print('\n *** Deformations will be calculated at: '+ctimeC+'\n')
    rdt = rTm2 - rTm1
    if not cfg.lc_accurate_time:
        print('      => time step to be used: `dt` = '+str(round(rdt,2))+' = '+str(round(rdt/cfg.rc_day2sec,2))+' days')

    vclck = split('_',ctimeC)
    chh, cmm = split(':',vclck[1])[0], split(':',vclck[1])[1]
    cclck = str.replace( vclck[0],'-','')+'-'+chh+'h'+cmm
    
    reskm, reskm2 = QUA1.reskm_nmnl, QUA2.reskm_nmnl
    if reskm != reskm2:
        print('ERROR: quads do not have the same nominal resolution in the 2 files:',reskm, reskm2)
        exit(0)
    print('\n *** Nominal resolution for the Quads in both files =',reskm,'km')
    
    # Comprehensive name for npz and figs to save later on:
    corigin = QUA1.origin    
    cfnm    = corigin

    # Some info from npz file name: #fixme: too much dependency on file name...
    cf1, cf2 = path.basename(cf_Q1), path.basename(cf_Q2)
    if corigin == 'RGPS':
        cbatch =  split('_',cf1)[2]
        cdtbin =  '_'+split('_',cf1)[3]
    elif split('_',corigin)[0] == 'NEMO-SI3':
        cbatch =  split('_',cf1)[-5]
        cdtbin = '_'+split('_',cf1)[-4]        
    else:
        print('FIXME: unknow origin: ',corigin); exit(0)
    if not cdtbin[1:3] in ['dt','No']:
        print('ERROR: we could not figure out `cdtbin`!'); exit

    #print('LOLO: cdtbin[1:3] =',cdtbin[1:3]);exit(0)
        
    # Test for coarsening realisations:
    cr1, cr2 = '', ''
    cc1, cc2 = split('-',split('_',cf1)[-1]), split('-',split('_',cf2)[-1])
    if len(cc1)==2: cr1 = cc1[0]+'-'
    if len(cc2)==2: cr2 = cc2[0]+'-'
    if cr1!='' and cr1!=cr2:
        print('ERROR: with the realisation coarsening res:',cr1,cr2); exit(0)
    else:
        print('     => realisation with `rd_ss` =',cc1[0],'km !!!')

    cfnm += '_'+cbatch+cdtbin
    cfnm += '_'+cclck        
    cfnm += '_'+cr1+str(reskm)+'km'
    #print('LOLO: cfnm =',cfnm); exit(0)
    
    print('\n *** Number of points in the two records:', QUA1.nP, QUA2.nP)
    print('\n *** Number of quads in the two records:' , QUA1.nQ, QUA2.nQ)

    if  cdtbin=='_NoBin':
        dtbin = 0
    else:
        dtbin = int(cdtbin[3:])*3600
        print('\n *** width of time bin used in RGPS =',dtbin/3600,'hours!')


    #kk = cfg.updateConfig4Scale( reskm )
    #print('\n *** Max time deviation accepted for vertices: `rc_t_dev_cancel` =',cfg.rc_t_dev_cancel,'s')    
    
    if cfg.lc_accurate_time:
        figSfx='_tbuoy.png'
    else:
        figSfx='_tglob.png'
    
    # The Quads we retain, i.e. those who exist in both snapshots:
    vnm, vidx1, vidx2 = np.intersect1d( QUA1.QuadNames, QUA2.QuadNames, assume_unique=True, return_indices=True )
    nQ = len(vnm) ; # also = len(vidx*)

    znm, zidx1, zidx2 = np.intersect1d( QUA1.QuadIDs, QUA2.QuadIDs, assume_unique=True, return_indices=True )
    nQ2 = len(znm)
    if nQ!=nQ2 or np.sum(zidx1-vidx1)!=0 or np.sum(zidx2-vidx2)!=0:
        print('ERROR: we do not get the same info based on Quad names and Quad Ids !!!')
        exit(0)

    print('       => there are '+str(nQ)+' Quads common to the 2 records!\n')


    # Now for some weird reasons time of a given buoy can be the same in the 2 quads:
    zTime1 = QUA1.MeshVrtcPntTime()[vidx1,:]
    zTime2 = QUA2.MeshVrtcPntTime()[vidx2,:]
    zdT = zTime2-zTime1
    #
    if np.any(zdT==0.):
        print('\n WARNING: time for some buoys is the same in the 2 records!')
        (idxKeep,) = np.where( (zdT[:,0]>0.) & (zdT[:,1]>0.) & (zdT[:,2]>0.) & (zdT[:,3]>0.) )
        idxKeep = np.array( idxKeep , dtype=int )
        nQn = len(idxKeep)
        print('   ==> '+str(nQ-nQn)+' quads /'+str(nQ)+' disregarded because they have the same time in both record !')
        vidx1 = vidx1[idxKeep]
        vidx2 = vidx2[idxKeep]
        nQ = nQn
    #
    del zTime1, zTime2, zdT


    # Coordinates of the 4 points of quadrangles for the 2 consecutive records:
    zXY1 = QUA1.MeshPointXY[vidx1,:,:].copy() ; #  km !
    zXY2 = QUA2.MeshPointXY[vidx2,:,:].copy() ; #  km !

    # Computation of partial derivative of velocity vector constructed from the 2 consecutive positions:
    if cfg.lc_accurate_time:
        # Time of poins of the 4 points of quadrangles for the 2 consecutive records:
        zTime1 = QUA1.MeshVrtcPntTime()[vidx1,:]
        zTime2 = QUA2.MeshVrtcPntTime()[vidx2,:]
        #
        zX, zY, zU, zV, zdUdxy, zdVdxy, zAq = mjt.PDVfromPos( 1, zXY1, zXY2, QUA1.area()[vidx1], QUA2.area()[vidx2], 
                                                              xtime1=zTime1, xtime2=zTime2, iverbose=idebug )
        # => having rdt=1 will yield fuck-up fields if used, and must not be used!
        
    else:

        zX, zY, zU, zV, zdUdxy, zdVdxy, zAq = mjt.PDVfromPos( rdt, zXY1, zXY2, QUA1.area()[vidx1], QUA2.area()[vidx2],
                                                              iverbose=idebug )

    # zX, zY => positions  of the 4 vertices at center of time interval!
    # zU, zV => velocities of the 4 vertices at center of time interval!
    #   zAq  => Quad area used in the estimation of `zdUdxy, zdVdxy` !

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
        mjt.ShowDefQuad( zXc, zYc, zUc, cfig=cdir+'/zvUc_'+cfnm+figSfx, cwhat='Uc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, unit='m/s' )
        mjt.ShowDefQuad( zXc, zYc, zVc, cfig=cdir+'/zvVc_'+cfnm+figSfx, cwhat='Vc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, unit='m/s' )
        mjt.ShowDefQuad( zXc, zYc, np.sqrt(zUc*zUc+zVc*zVc), cfig=cdir+'/zUMc_'+cfnm+figSfx, cwhat='UMc',
                             pFmin=0., pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, unit='m/s' )

        #mjt.ShowDefQuad( zX, zY, cfg.rc_day2sec*zdiv, cfig=cdir+'/zd_'+cfnm+'_Divergence'+figSfx, cwhat='div',
        #                 pFmin=-cfg.rc_div_max, pFmax=cfg.rc_div_max, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
        #                 title=corigin+': divergence' )


        
        #
        del zUc, zVc

    # Divergence aka SigmaI:
    zdiv = mjt.DivPDV(zdUdxy, zdVdxy)

    # `Maximum Shear Strain Rate`**2:
    zshr2 = mjt.Shear2PDV(zdUdxy, zdVdxy )

    del zdUdxy, zdVdxy

    # Maximum Shear Strain Rate aka SigmaII:
    zshr = np.zeros(np.shape(zshr2))
    zshr = np.sqrt(zshr2)

    # Total deformation rate:
    ztot = np.zeros(np.shape(zdiv))
    ztot[:] = np.sqrt( zdiv[:]*zdiv[:] + zshr2[:] )

    del zshr2


    # Saving data:
    np.savez_compressed( './npz/DEFORMATIONS_'+cfnm+'.npz', time=rtimeC, date=ctimeC, Npoints=nQ,
                         Xc=zXc, Yc=zYc, divergence=zdiv, shear=zshr, total=ztot, quadArea=zAq, origin=corigin, reskm_nmnl=reskm )


    # Some plots:
    if iplot>0:

        if corigin != 'RGPS':
            corigin = str.replace( corigin,'NEMO-','')
            corigin = str.replace( corigin,'_NANUK4_','-')
        
        cresinfo = '('+str(reskm)+' km)'

        nmproj=NameArcticProj
        # Filled quads projected on the Arctic map:
        mjt.ShowDefQuadGeoArctic( zX, zY, cfg.rc_day2sec*zdiv, cfig=cdir+'/map_zd_'+cfnm+'_Divergence'+figSfx, nmproj=NameArcticProj, cwhat='div',
                                  pFmin=-cfg.rc_div_max, pFmax=cfg.rc_div_max, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': divergence '+cresinfo )

        mjt.ShowDefQuadGeoArctic( zX, zY, cfg.rc_day2sec*zshr, cfig=cdir+'/map_zs_'+cfnm+'_Shear'+figSfx,      nmproj=NameArcticProj, cwhat='shr',
                                  pFmin=0.,      pFmax=cfg.rc_shr_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': shear '+cresinfo )

        mjt.ShowDefQuadGeoArctic( zX, zY, cfg.rc_day2sec*zshr, cfig=cdir+'/map_zt_'+cfnm+'_Total'+figSfx,      nmproj=NameArcticProj, cwhat='tot',
                                  pFmin=0.,      pFmax=cfg.rc_tot_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': total deformation '+cresinfo )

    if iplot>1:
        # Filled quads projected on RGPS projection (Cartesian):
        mjt.ShowDefQuad( zX, zY, cfg.rc_day2sec*zdiv, cfig=cdir+'/zd_'+cfnm+'_Divergence'+figSfx, cwhat='div',
                                  pFmin=-cfg.rc_div_max, pFmax=cfg.rc_div_max, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': divergence '+cresinfo )

        mjt.ShowDefQuad( zX, zY, cfg.rc_day2sec*zshr, cfig=cdir+'/zs_'+cfnm+'_Shear'+figSfx,      cwhat='shr',
                                  pFmin=0.,      pFmax=cfg.rc_shr_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': shear '+cresinfo )

        mjt.ShowDefQuad( zX, zY, cfg.rc_day2sec*zshr, cfig=cdir+'/zt_'+cfnm+'_Total'+figSfx,      cwhat='tot',
                                  pFmin=0.,      pFmax=cfg.rc_tot_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': total deformation '+cresinfo )


        ###
