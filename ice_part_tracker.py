#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

#TODO:


from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

import gonzag as gz

from climporn import epoch2clock, clock2epoch
import mojito   as mjt



idebug=2




if __name__ == '__main__':

    if not len(argv) in [4,5]:
        print('Usage: '+argv[0]+' <ice_file_velocities_SI3> <mesh_mask>')
        exit(0)
    cf_Q1 = argv[1]
    cf_Q2 = argv[2]
    time_dev_max= float(argv[3])
    #
    marker_size=None
    if len(argv) == 5:
        marker_size = int(argv[4])

    # Reading the quad meshes in both npz files:
    QUA1 = mjt.LoadClassPolygon( cf_Q1, ctype='Q' )
    QUA2 = mjt.LoadClassPolygon( cf_Q2, ctype='Q' )

    # Debug have a look at the times of all points and get the actual mean time for each file:
    rT1 = CheckTimeSanityQuad(1, QUA1, time_dev_max, iverbose=idebug)
    rT2 = CheckTimeSanityQuad(2, QUA2, time_dev_max, iverbose=idebug)

    rtimeC = 0.5*(rT1+rT2)
    ctimeC = epoch2clock(rtimeC)
    print('\n *** Deformations will be calculated at: '+ctimeC+'\n')
    rdt = rT2 - rT1
    print('      => time step to be used: `dt` = '+str(round(rdt,2))+' = '+str(round(rdt/(3600*24),2))+' days')

    vclck = split('_',ctimeC)
    chh = split(':',vclck[1])[0]
    cclck = str.replace( vclck[0],'-','')+'-'+chh+'h'


    # Comprehensive name for npz and figs to save later on:
    cf1, cf2 = path.basename(cf_Q1), path.basename(cf_Q2)
    cfnm = split('_',cf1)[1]
    cfnm = cfnm+'_'+cclck

    # Try to get a spatial resolution scale from the name:
    cres = ''
    cc1 = split('km_',cf1)[0] ; cc1 = split('_',cc1)[-1]
    cc2 = split('km_',cf2)[0] ; cc2 = split('_',cc2)[-1]
    if cc1==cc2 and cc1.isdigit() and cc2.isdigit():
        cres = '_'+cc1+'km'

    print('\n *** Number of points in the two records:',QUA1.nP,QUA2.nP)
    print('\n *** Number of quads in the two records:',QUA1.nQ,QUA2.nQ)

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



    if idebug>0:
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
        mjt.ShowDeformation( zzx, zzy, zzu, cfig='./figs/zvU4_'+cfnm+'_'+cres+'.png', cwhat='U4',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zzx, zzy, zzv, cfig='./figs/zvV4_'+cfnm+'_'+cres+'.png', cwhat='V4',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zzx, zzy, np.sqrt(zzu*zzu+zzv*zzv), cfig='./figs/zUM4_'+cfnm+'_'+cres+'.png', cwhat='UMc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        #
        del zzx, zzy, zzu, zzv

    # Coordinates of barycenter of Quads at center of time interval:
    zXc = np.mean( zX[:,:], axis=1 )
    zYc = np.mean( zY[:,:], axis=1 )

    if idebug>1:
        # Velocities of barycenter of Quads at center of time interval:
        zUc = np.mean( zU[:,:], axis=1 )
        zVc = np.mean( zV[:,:], axis=1 )
        #
        mjt.ShowDeformation( zXc, zYc, zUc, cfig='./figs/zvUc_'+cfnm+'_'+cres+'.png', cwhat='Uc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zXc, zYc, zVc, cfig='./figs/zvVc_'+cfnm+'_'+cres+'.png', cwhat='Vc',
                             pFmin=-0.2, pFmax=0.2, zoom=zoom, rangeX=zrx, rangeY=zry, marker_size=marker_size, unit='m/s' )
        mjt.ShowDeformation( zXc, zYc, np.sqrt(zUc*zUc+zVc*zVc), cfig='./figs/zUMc_'+cfnm+'_'+cres+'.png', cwhat='UMc',
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
                         Xc=zXc, Yc=zYc, divergence=zdiv, shear=zshr, total=ztot )


    # Some plots:
    if not path.exists('./figs'): mkdir('./figs')

    mjt.ShowDeformation( zXc, zYc, rconv*zdiv, cfig='./figs/zd_'+cfnm+'_Divergence'+cres+'.png', cwhat='div',
                         pFmin=-div_max, pFmax=div_max, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$', marker_size=marker_size )
    mjt.ShowDeformation( zXc, zYc, rconv*zshr, cfig='./figs/zs_'+cfnm+'_Shear'+cres+'.png',      cwhat='shr',
                         pFmin=0.,      pFmax=shr_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$', marker_size=marker_size )
    mjt.ShowDeformation( zXc, zYc, rconv*zshr, cfig='./figs/zt_'+cfnm+'_Total'+cres+'.png',      cwhat='tot',
                         pFmin=0.,      pFmax=tot_max,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$', marker_size=marker_size )


    ###
