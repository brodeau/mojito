#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
from glob import glob
import numpy as np
from re import split
#
import mojito   as mjt
from mojito import config as cfg

idebug=1
iplot=1

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

if __name__ == '__main__':

    kk = cfg.initialize()
    
    if not len(argv) in [6]:
        print('Usage: '+argv[0]+' <directory_input_npz_files> <batch.N> <dtbin_h> <creskm> <string_id_origin>')
        exit(0)    

    cd_in  = argv[1]
    kbatch = argv[2]
    cdtbin = argv[3]
    creskm = argv[4]
    cidorg = argv[5]
    
    sbatch = 'S%3.3i'%(int(kbatch))

    
    
    print(cd_in+'/'+cprefixIn+'*'+cidorg+'_'+sbatch+'_'+cdtbin+'*_'+creskm+'km.npz')
    listnpz1 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'_'+sbatch+'_'+cdtbin+'*_'+creskm+'km.npz') )
    listnpz2 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'_'+sbatch+'_'+cdtbin+'*_*-'+creskm+'km.npz') )

    if len(listnpz1)>0 and len(listnpz2)>0:
        print('ERROR: we have both npz files with suffixes lile `*_Xkm` and `_Y-Xkm` !!!'); exit(0)
    if len(listnpz2)>0:
        lrlstKM = True
        listnpz = listnpz2
    else:
        lrlstKM = False
        listnpz = listnpz1
        
    
    # Polpulating deformation files available:    
    nbF = len(listnpz)
    print('\n *** We found '+str(nbF)+' deformation files into '+cd_in+' !')
    print('      => files to gather into a single one:',listnpz,'\n')

    # First, populate the number of deformation quadrangles in each file:
    kNbPoints = np.zeros(nbF, dtype=int ) ; # number of points in file
    vsubRes   = np.zeros(nbF, dtype=int ) ; # coarsifying resolution

    
    kf = 0
    for ff in listnpz:

        css   = split( '-', split( '_', split('\.',ff)[0] )[-1] )[0]
        vsubRes[kf] = int(css)
        csskm = css+'km'
        
        with np.load(ff) as data:
            nPnts =      data['Npoints']
        print('\n  # File: '+ff+' ('+csskm+') => '+str(nPnts)+' deformation quads!')
        kNbPoints[kf] = nPnts
        kf+=1

    
    # Winner:
    kw = np.argmax(kNbPoints)
    Nmax = kNbPoints[kw]
    (idxmax,) = np.where(kNbPoints==Nmax)
    nn = len(idxmax)    
    # Takes the one the closest to the reference resolution
    if nn>1:
        ki = np.argmin( np.abs(vsubRes[idxmax] - int(creskm)) )
        print(vsubRes[idxmax] - (int(creskm)-2), ki) ; # "-2" to be closer to the actual scale used to coarsify!
        kw = idxmax[ki]
    Nmax = kNbPoints[kw]
    print(' *** Winner is at index:',kw,'(',Nmax,'points),',vsubRes[kw],'km,',nn,'files with this number of points.')



    # The "in use" quads coordinates array:
    nQmaxF = nbF*Nmax ; # the ever max number of quads we expect in a single file, final version of `nQmaxF` will have less points...    
    zInUseXY = np.zeros((nQmaxF,4,2)) - 9999.
    zmsk     = np.zeros( nQmaxF , dtype='i1' )
    
    # 1/ Save all the quads of winner:
    with np.load(listnpz[kw]) as data:
        #    zInUseXY = data[
        #zX4 = data['X4']
        #zY4 = data['Y4']
        #print('shape zX4 =',np.shape(zX4))
        zInUseXY[0:Nmax,:,0], zInUseXY[0:Nmax,:,1] = data['X4'], data['Y4']

    #print(zX4)
        
    

        
    #zX4 = np.array(zX4)
    #print(np.shape(zX4))
    
    exit(0)


    
    kBatchName  = np.zeros(nbF, dtype='U4')
    kiDate      = np.zeros(nbF, dtype=int ) ; # date in epoch time at which deformations were calculated
    kNbPoints   = np.zeros(nbF, dtype=int ) ; # number of points in file    
    
    list_date = []
    kf = 0
    for ff in listnpz:
        print('\n  # File: '+ff)

        with np.load(ff) as data:
            rdate = int( data['time'] )
            nPnts =      data['Npoints']
            if kf==0:
                corigin = str(data['origin'])
                reskm   = int(data['reskm_nmnl'])
        fb = path.basename(ff)
        vf = split('_',fb)


        cdateh = split('-',vf[-2])[0]
        list_date.append(split('h',cdateh)[0])
        
        kBatchName[kf] = vf[-4]
    
        kiDate[kf] = rdate
        kNbPoints[kf] = nPnts
    
        print('   * Batch: '+kBatchName[kf] )
        print('   * Date = ',mjt.epoch2clock(kiDate[kf]))
        print('   * Nb. of points = ',kNbPoints[kf] )

        kf = kf+1
    
    print('\n')

    nP = np.sum(kNbPoints)
    print('  ==> Total number of points:', nP)
    print('  ==> list of dates:', list_date[:])
    
    cdt1, cdt2 = list_date[0],list_date[-1]
    cperiod = cdt1+'-'+cdt2

    dtbin=0
    cinfobin = ''
    if cdtbin!='nemoTsi3_NoBin':
        dtbin=int(cdtbin)
        cinfobin = '_dt'+cdtbin
        
    if str(reskm) != creskm:
        print('ERROR: spatial scale (km) passed as argument does not match that found in deformation files!')
        exit(0)
    
    # Now that we know the total number of points we can allocate and fill arrays for divergence and shear    
    Zshr = np.zeros(nP)
    ZDiv = np.zeros(nP)
    Ztot = np.zeros(nP)
    ZAqd = np.zeros(nP)
    Zdat = np.zeros(nP,dtype=int)
    
    jP = 0 ; # Counter from 0 to nP-1
    kf = 0 ; # Counter for files, 0 to nbF-1
    for ff in listnpz:
        jPe = jP + kNbPoints[kf]
        with np.load(ff) as data:
            zdiv  =      data['divergence']
            zshr  =      data['shear']
            ztot  =      data['total']
            za    =      data['quadArea']
        #kf
        ZDiv[jP:jPe] = cfg.rc_day2sec*zdiv ; # day^-1
        Zshr[jP:jPe] = cfg.rc_day2sec*zshr ; # day^-1
        Ztot[jP:jPe] = cfg.rc_day2sec*ztot ; # day^-1
        ZAqd[jP:jPe] =     za     ; # km^2
        Zdat[jP:jPe] =  kiDate[kf]
        #
        jP = jPe
        kf = kf+1

    del zdiv, zshr, ztot, za

    
    cfroot = corigin+cinfobin+'_'+str(reskm)+'km_'+cperiod   

    
    # Save it for the scaling and PDFs to come ...
    print(' *** Saving: '+cd_in+'/def_SHR_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_SHR_'+cfroot+'.npz', name='shear', origin=corigin,
                         Nbatch=nbF, dates_batch=kiDate,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nP, xshr=Zshr, quadArea=ZAqd, dates_point=Zdat )

    print(' *** Saving: '+cd_in+'/def_DIV_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_DIV_'+cfroot+'.npz', name='divergence', origin=corigin,
                         Nbatch=nbF, dates_batch=kiDate,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nP, xdiv=ZDiv, quadArea=ZAqd, dates_point=Zdat )

    print(' *** Saving: '+cd_in+'/def_TOT_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_TOT_'+cfroot+'.npz', name='shear', origin=corigin,
                         Nbatch=nbF, dates_batch=kiDate,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nP, xtot=Ztot, quadArea=ZAqd, dates_point=Zdat )

