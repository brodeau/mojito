#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

#TODO:
#     * figure out what the acceptable "time_dev_from_mean_allowed" and if some points go beyond it, just remove them from the data!


from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=2
iplot=1

l_accurate_time=True

# Conversion from s-1 to day-1:
rconv = 24.*3600.

# For figures, in days^-1:
div_max = 0.1
shr_max = 0.1
tot_max = 0.1

zoom=1


if __name__ == '__main__':

    if not len(argv) in [4,5]:
        print('Usage: '+argv[0]+' <file_Q_mesh_N1.npz> <file_Q_mesh_N2.npz> <time_dev_from_mean_allowed (s)> (<marker_size>)')
        exit(0)
    cf_Q1 = argv[1]
    cf_Q2 = argv[2]
    time_dev_max= float(argv[3])
    #
    marker_size=None
    if len(argv) == 5:
        marker_size = int(argv[4])


    if iplot>0:
        cdir = './figs/deformation'
        if not path.exists(cdir): mkdir(cdir)
        
    # Reading the quad meshes in both npz files:
    QUA1 = mjt.LoadClassPolygon( cf_Q1, ctype='Q' )
    QUA2 = mjt.LoadClassPolygon( cf_Q2, ctype='Q' )

    # Debug have a look at the times of all points and get the actual mean time for each file:
    rT1 = mjt.CheckTimeConsistencyQuads(1, QUA1, time_dev_max, iverbose=idebug)
    rT2 = mjt.CheckTimeConsistencyQuads(2, QUA2, time_dev_max, iverbose=idebug)

    rtimeC = 0.5*(rT1+rT2)
    ctimeC = epoch2clock(rtimeC)
    print('\n *** Deformations will be calculated at: '+ctimeC+'\n')
    rdt = rT2 - rT1
    print('      => time step to be used: `dt` = '+str(round(rdt,2))+' = '+str(round(rdt/(3600*24),2))+' days')

    vclck = split('_',ctimeC)
    chh = split(':',vclck[1])[0]
    cclck = str.replace( vclck[0],'-','')+'h'+chh

    
    # Comprehensive name for npz and figs to save later on:
    corigin = QUA1.origin
    cfnm = corigin+'_'+cclck

    creskm = ''
        
    # Some info from npz file name:
    cf1, cf2 = path.basename(cf_Q1), path.basename(cf_Q2)
    if corigin == 'RGPS':
        cc1 = split('km.',cf1)[0] ; cc1 = split('_',cc1)[-1]
        cc2 = split('km.',cf2)[0] ; cc2 = split('_',cc2)[-1]
    else:    
        cc1 = split('km_',cf1)[0] ; cc1 = split('_',cc1)[-1]
        cc2 = split('km_',cf2)[0] ; cc2 = split('_',cc2)[-1]        
    if cc1==cc2 and cc1.isdigit() and cc2.isdigit():
        creskm = '_'+cc1+'km'
        print('\n *** We have a resolution from names of files => '+cc1+' km !')
        cfnm += creskm
        
    print('\n *** Number of points in the two records:', QUA1.nP, QUA2.nP)
    print('\n *** Number of quads in the two records:' , QUA1.nQ, QUA2.nQ)
    
    # The Quads we retain, i.e. those who exist in both snapshots:
    vnm, vidx1, vidx2 = np.intersect1d( QUA1.QuadNames, QUA2.QuadNames, assume_unique=True, return_indices=True )
    nQ = len(vnm) ; # also = len(vidx*)

    znm, zidx1, zidx2 = np.intersect1d( QUA1.QuadIDs, QUA2.QuadIDs, assume_unique=True, return_indices=True )
    nQ2 = len(znm)
    if nQ!=nQ2 or np.sum(zidx1-vidx1)!=0 or np.sum(zidx2-vidx2)!=0:
        print('ERROR: we do not get the same info based on Quad names and Quad Ids !!!')
        exit(0)

    print('       => there are '+str(nQ)+' Quads common to the 2 records!\n')

    # Coordinates of the 4 points of quadrangles for the 2 consecutive records:
    zXY1 = QUA1.MeshPointXY[vidx1,:,:].copy() ; #  km !
    zXY2 = QUA2.MeshPointXY[vidx2,:,:].copy() ; #  km !

    # Computation of partial derivative of velocity vector constructed from the 2 consecutive positions:
    if l_accurate_time:
        # Time of poins of the 4 points of quadrangles for the 2 consecutive records:
        zTime1 = QUA1.MeshVrtcPntTime()[vidx1,:]
        zTime2 = QUA2.MeshVrtcPntTime()[vidx2,:]
        #
        zX, zY, zU, zV, zdUdxy, zdVdxy = mjt.PDVfromPos( rdt, zXY1, zXY2, QUA1.area()[vidx1], QUA2.area()[vidx2],
                                                         xtime1=zTime1, xtime2=zTime2, iverbose=idebug )
    else:

        zX, zY, zU, zV, zdUdxy, zdVdxy = mjt.PDVfromPos( rdt, zXY1, zXY2, QUA1.area()[vidx1], QUA2.area()[vidx2],
                                                         iverbose=idebug )

    # zX, zY => positions  of the 4 vertices at center of time interval!
    # zU, zV => velocities of the 4 vertices at center of time interval!

    zrx = [ np.min(zX)-25. , np.max(zX)+25. ]
    zry = [ np.min(zY)-25. , np.max(zY)+25. ]



    if idebug>0 and iplot>0:
        # DEBUG: show velocities at each 4 vertices of each Quad:
        zzx = zX.flatten()
        zzy = zY.flatten()
        zzu = zU.flatten()
        zzv = zV.flatten()
        #
        zcoor0 = np.array([ zzx, zzy ]).T        
        _,idx_uniq = np.unique(zcoor0, axis=0, return_index=True) ; # index location to get rid of clones... (because lot of the same points are used by the same quads)
        del zcoor0

        zzx = zzx[idx_uniq]
        zzy = zzy[idx_uniq]
        zzu = zzu[idx_uniq]
        zzv = zzv[idx_uniq]
        #
        mjt.ShowDeformation( zzx, zzy, zzu, cfig=cdir+'/zvU4_'+cfnm+'.png', cwhat='U4',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zzx, zzy, zzv, cfig=cdir+'/zvV4_'+cfnm+'.png', cwhat='V4',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zzx, zzy, np.sqrt(zzu*zzu+zzv*zzv), cfig=cdir+'/zUM4_'+cfnm+'.png', cwhat='UMc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        #
        del zzx, zzy, zzu, zzv

    # Coordinates of barycenter of Quads at center of time interval:
    zXc = np.mean( zX[:,:], axis=1 )
    zYc = np.mean( zY[:,:], axis=1 )

    if idebug>1 and iplot>0:
        # Velocities of barycenter of Quads at center of time interval:
        zUc = np.mean( zU[:,:], axis=1 )
        zVc = np.mean( zV[:,:], axis=1 )
        #
        mjt.ShowDeformation( zXc, zYc, zUc, cfig=cdir+'/zvUc_'+cfnm+'.png', cwhat='Uc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zXc, zYc, zVc, cfig=cdir+'/zvVc_'+cfnm+'.png', cwhat='Vc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zXc, zYc, np.sqrt(zUc*zUc+zVc*zVc), cfig=cdir+'/zUMc_'+cfnm+'.png', cwhat='UMc',
                             pFmin=0., pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
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
                         Xc=zXc, Yc=zYc, divergence=zdiv, shear=zshr, total=ztot, origin=corigin )


    # Some plots:
    if iplot>0:
        mjt.ShowDeformation( zXc, zYc, rconv*zdiv, cfig=cdir+'/zd_'+cfnm+'_Divergence.png', cwhat='div',
                             pFmin=-div_max, pFmax=div_max, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                             marker_size=marker_size, title=corigin+': divergence' )
        mjt.ShowDeformation( zXc, zYc, rconv*zshr, cfig=cdir+'/zs_'+cfnm+'_Shear.png',      cwhat='shr',
                             pFmin=0.,      pFmax=shr_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                             marker_size=marker_size, title=corigin+': shear' )
        mjt.ShowDeformation( zXc, zYc, rconv*zshr, cfig=cdir+'/zt_'+cfnm+'_Total.png',      cwhat='tot',
                             pFmin=0.,      pFmax=tot_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                             marker_size=marker_size, title=corigin+': total deformation' )


    ###
