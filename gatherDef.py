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

    cslct = ''
    
    if not len(argv) in [5,6]:
        print('Usage: '+argv[0]+' <directory_input_npz_files> <dtbin_h> <creskm> <string_id_origin> (<S>)')
        exit(0)    

    cd_in  = argv[1]
    cdtbin = argv[2]
    creskm = argv[3]
    cidorg = argv[4]

    if len(argv)==6:
        lgroup = (argv[5]=='S')
        if lgroup:
            listnpz = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'*_dt'+cdtbin+'*_SLCT'+creskm+'km.npz') )
            print(cd_in+'/'+cprefixIn+'*'+cidorg+'*_dt'+cdtbin+'*_SLCT'+creskm+'km.npz')
        else:
            mjt.printEE('argument 6 can only be "S"!!!')

    else:
        listnpz1 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'*'+cdtbin+'*_'+creskm+'km.npz') )
        listnpz2 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'*'+cdtbin+'*_*-'+creskm+'km.npz') )
    
        if len(listnpz1)>0 and len(listnpz2)>0:
            mjt.printEE('we have both npz files with suffixes lile `*_Xkm` and `_Y-Xkm` !!!')
        if len(listnpz2)>0:
            lrlstKM = True
            listnpz = listnpz2
        else:
            lrlstKM = False
            listnpz = listnpz1
        
                    
    # Polpulating deformation files available:    
    nbFiles = len(listnpz)
    print('\n *** We found '+str(nbFiles)+' deformation files into '+cd_in+' !')
    print('      => files to gather into a single one:',listnpz,'\n')

    kBatchName  = np.zeros(nbFiles, dtype='U4')
    #kiDate      = np.zeros(nbFiles, dtype=int ) ; # date in epoch time at which deformations were calculated
    kNbPoints   = np.zeros(nbFiles, dtype=int ) ; # number of points in file    
    
    #list_date = []
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
        #cdateh = split('-',vf[-2])[0]
        #list_date.append(split('h',cdateh)[0])

        if lgroup:
            kBatchName[kf] = vf[-3]
        else:
            kBatchName[kf] = vf[-4]
    
        #kiDate[kf] = rdate
        kNbPoints[kf] = nPnts
    
        print('   * Batch: '+kBatchName[kf] )
        #print('   * Date = ',mjt.epoch2clock(kiDate[kf]))
        print('   * Nb. of points = ',kNbPoints[kf] )

        kf = kf+1
    
    print('\n')
    #exit(0)
    nP = np.sum(kNbPoints)
    print('  ==> Total number of points:', nP)
    #print('  ==> list of dates:', list_date[:])
    #cdt1, cdt2 = list_date[0],list_date[-1]
    #cperiod = cdt1+'-'+cdt2

    dtbin=0
    cinfobin = ''
    if cdtbin!='nemoTsi3_NoBin':
        dtbin=int(cdtbin)
        cinfobin = '_dt'+cdtbin
        
    if str(reskm) != creskm:
        mjt.printEE('spatial scale (km) passed as argument does not match that found in deformation files!')
    
    # Now that we know the total number of points we can allocate and fill arrays for divergence and shear    
    Zshr = np.zeros(nP)
    ZDiv = np.zeros(nP)
    Ztot = np.zeros(nP)
    ZAqd = np.zeros(nP)
    #Zdat = np.zeros(nP,dtype=int)
    
    jP = 0 ; # Counter from 0 to nP-1
    kf = 0 ; # Counter for files, 0 to nbFiles-1
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
        #Zdat[jP:jPe] =  kiDate[kf]
        #
        jP = jPe
        kf = kf+1

    del zdiv, zshr, ztot, za

    
    #cfroot = corigin+cinfobin+'_'+str(reskm)+'km_'+cperiod
    cfroot = corigin+cinfobin+'_'+str(reskm)+'km'

    
    # Save it for the scaling and PDFs to come ...
    #   , dates_batch=kiDate, period=cperiod
    print(' *** Saving: '+cd_in+'/def_SHR_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_SHR_'+cfroot+'.npz', name='shear', origin=corigin,
                         Nbatch=nbFiles,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nP, xshr=Zshr, quadArea=ZAqd) ; #, dates_point=Zdat )

    print(' *** Saving: '+cd_in+'/def_DIV_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_DIV_'+cfroot+'.npz', name='divergence', origin=corigin,
                         Nbatch=nbFiles,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nP, xdiv=ZDiv, quadArea=ZAqd) ; #, dates_point=Zdat )

    print(' *** Saving: '+cd_in+'/def_TOT_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_TOT_'+cfroot+'.npz', name='shear', origin=corigin,
                         Nbatch=nbFiles,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nP, xtot=Ztot, quadArea=ZAqd) ; #, dates_point=Zdat )

