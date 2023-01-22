#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

#TODO:
#     * figure out what the acceptable "time_dev_from_mean_allowed" and if some points go beyond it, just remove them from the data!


from sys import argv, exit
from os import path
import numpy as np
from re import split

from scipy.spatial import Delaunay

from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=0

l_accurate_time=True


def CheckTimeSanityQuad( kF, QD, time_dev_from_mean_allowed ):
    '''
        * kF: file number
        * QD: Quad class loaded from file `kF`
    '''
    print('\n *** In file #'+str(kF)+':')
    print('     (time_dev_from_mean_allowed =', time_dev_from_mean_allowed/60.,' minutes)' )
    rTmean = np.mean(QD.PointTime)
    cTmean = epoch2clock(rTmean)
    rStdDv = mjt.StdDev(rTmean, QD.PointTime)
    print('     => Actual Mean time =',cTmean,', Standard Deviation =',rStdDv/60.,' minutes!')
    zadiff = np.abs(QD.PointTime-rTmean)
    zdt = np.max(zadiff)/60.
    print('     ==> furthest point is '+str(round(zdt,2))+' minutes away from mean!')

    if np.any(zadiff>time_dev_from_mean_allowed):
        print('ERROR: some points read in file #'+str(kF)+' are too far (in time) from mean time=',cTmean)
        exit(0)
    else:
        print('     => ok! No points further than '+str(time_dev_from_mean_allowed/60.)+' minutes from mean...')
        
    return rTmean



if __name__ == '__main__':

    if not len(argv) in [5]:
        print('Usage: '+argv[0]+' <file_Q_mesh_N1.npz> <file_Q_mesh_N2.npz> <time_dev_from_mean_allowed (s)> <marker_size>')
        exit(0)
    cf_Q1 = argv[1]
    cf_Q2 = argv[2]
    time_dev_max= float(argv[3])
    mrkrsz= int(argv[4])
        
    # Reading the quad meshes in both npz files:
    QUA1 = mjt.LoadClassPolygon( cf_Q1, ctype='Q' )
    QUA2 = mjt.LoadClassPolygon( cf_Q2, ctype='Q' )

    #cdt1 = str(QUA1.date)
    #cdt2 = str(QUA2.date)
    #print('\n *** Mean date for the two records:',cdt1,cdt2)
    #idt1 = clock2epoch(cdt1)
    #idt2 = clock2epoch(cdt2)
        
    # Debug have a look at the times of all points and get the actual mean time for each file:
    rT1 = CheckTimeSanityQuad(1, QUA1, time_dev_max)
    rT2 = CheckTimeSanityQuad(2, QUA2, time_dev_max)

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
    cnm_pref = split('_',cf1)[1]
    #ccdt1, ccdt2 = split('_',cdt1)[0], split('_',cdt2)[0]
    #cnm_pref = 'z_'+cnm_pref+'_'+ccdt1+'_'+ccdt2
    cnm_pref = 'z_'+cnm_pref+'_'+cclck
    
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
    zXY1 = QUA1.MeshPointXY[vidx1,:,:].copy() ; #* 1000.  ; # 1000 => from km to m
    zXY2 = QUA2.MeshPointXY[vidx2,:,:].copy() ; #* 1000.

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
    
    # zX, zY => positions at center of time interval!
    # zU, zV => velocities at center of time interval!
    
    
    zrx = [ np.min(zX)-25. , np.max(zX)+25. ]
    zry = [ np.min(zY)-25. , np.max(zY)+25. ]
    
    
    
    if idebug>0:
        # DEBUG: show velocities at each 4 vertices of each Quad:
        zzx = zX.flatten()
        zzy = zY.flatten()
        zzu = zU.flatten()
        zzv = zV.flatten()
        #
        mjt.ShowDeformation( zzx, zzy, zzu, cfig='./figs/'+cnm_pref+'_U4'+cres+'.png', cwhat='U4',
                             marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4, rangeX=zrx, rangeY=zry )
        mjt.ShowDeformation( zzx, zzy, zzv, cfig='./figs/'+cnm_pref+'_V4'+cres+'.png', cwhat='V4',
                             marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4, rangeX=zrx, rangeY=zry )
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
        mjt.ShowDeformation( zXc, zYc, zUc, cfig='./figs/'+cnm_pref+'_Uc'+cres+'.png', cwhat='Uc',
                             marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4, rangeX=zrx, rangeY=zry )
        mjt.ShowDeformation( zXc, zYc, zVc, cfig='./figs/'+cnm_pref+'_Vc'+cres+'.png', cwhat='Vc',
                             marker_size=mrkrsz, pFmin=-1e-4, pFmax=1.e-4, zoom=4, rangeX=zrx, rangeY=zry )
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
    np.savez_compressed( './npz/DEFORMATIONS_'+cnm_pref+'.npz', time=rtimeC, date=ctimeC, Npoints=nQ,
                         Xc=zXc, Yc=zYc, divergence=zdiv, shear=zshr, total=ztot )
    
    
    # Some plots:
    if not path.exists('./figs'): mkdir('./figs')
    
    mjt.ShowDeformation( zXc, zYc, zdiv, cfig='./figs/'+cnm_pref+'_Divergence'+cres+'.png', cwhat='div',
                         marker_size=mrkrsz, pFmin=-1.e-6, pFmax=1.e-6, zoom=4, rangeX=zrx, rangeY=zry )
    mjt.ShowDeformation( zXc, zYc, zshr, cfig='./figs/'+cnm_pref+'_Shear'+cres+'.png',      cwhat='shr',
                         marker_size=mrkrsz, pFmin=0.,      pFmax=0.8e-6,  zoom=4, rangeX=zrx, rangeY=zry )
    mjt.ShowDeformation( zXc, zYc, zshr, cfig='./figs/'+cnm_pref+'_Total'+cres+'.png',      cwhat='tot',
                         marker_size=mrkrsz, pFmin=0.,      pFmax=0.8e-6,  zoom=4, rangeX=zrx, rangeY=zry )
    
    
    ###
