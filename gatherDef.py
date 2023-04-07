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

idebug=1
iplot=1

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

# Conversion from s-1 to day-1:
rconv = 24.*3600.


if __name__ == '__main__':
    
    if not len(argv) in [5]:
        print('Usage: '+argv[0]+' <directory_input_npz_files> <dtbin_h> <creskm> <string_id_origin>')
        exit(0)    

    cd_in  = argv[1]
    cdtbin = argv[2]
    creskm = argv[3]
    cidorg = argv[4]

    
    listnpz1 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'*_dt'+cdtbin+'*_'+creskm+'km.npz') )
    listnpz2 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'*_dt'+cdtbin+'*_*-'+creskm+'km.npz') )

    if len(listnpz1)>0 and len(listnpz2)>0:
        print('ERROR: we have both npz files with suffixes lile `*_Xkm` and `_Y-Xkm` !!!'); exit(0)
    if len(listnpz2)>0:
        lrlstKM = True
        listnpz = listnpz2
    else:
        lrlstKM = False
        listnpz = listnpz1
        
    
    #exit(0)
    
    # Polpulating deformation files available:    
    nbFiles = len(listnpz)
    print('\n *** We found '+str(nbFiles)+' deformation files into '+cd_in+' !')
    print('      => files to gather into a single one:',listnpz,'\n')
    
    kBatchName = np.zeros(nbFiles, dtype='U4')
    kiDate      = np.zeros(nbFiles, dtype=int ) ; # date in epoch time at which deformations were calculated
    kNbPoints   = np.zeros(nbFiles, dtype=int ) ; # number of points in file    
    
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

    dtbin=int(cdtbin)
    if str(reskm) != creskm:
        print('ERROR: spatial scale (km) passed as argument does not match that found in deformation files!')
        exit(0)
    
    # Now that we know the total number of points we can allocate and fill arrays for divergence and shear    
    Zshr = np.zeros(nP)
    ZDiv = np.zeros(nP)
    Ztot = np.zeros(nP)
    ZAqd = np.zeros(nP)
    
    jP = 0 ; # Counter from 0 to nP-1
    kf = 0 ; # Counter for files, 0 to nbFiles-1
    for ff in listnpz:
        jPe = jP + kNbPoints[kf]
        with np.load(ff) as data:
            zdiv  =      data['divergence']
            zshr  =      data['shear']
            ztot  =      data['total']
            za    =      data['quadArea']
        #
        ZDiv[jP:jPe] = rconv*zdiv ; # day^-1
        Zshr[jP:jPe] = rconv*zshr ; # day^-1
        Ztot[jP:jPe] = rconv*ztot ; # day^-1
        ZAqd[jP:jPe] =     za     ; # km^2
        #
        jP = jPe
        kf = kf+1

    del zdiv, zshr, ztot, za

    
    cfroot = corigin+'_dt'+cdtbin+'_'+str(reskm)+'km_'+cperiod   

    
    # Save it for the scaling and PDFs to come ...
    print(' *** Saving: '+cd_in+'/def_SHR_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_SHR_'+cfroot+'.npz', name='shear', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin, Nbatch=nbFiles,
                         Np=nP, xshr=Zshr, quadArea=ZAqd )

    print(' *** Saving: '+cd_in+'/def_DIV_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_DIV_'+cfroot+'.npz', name='divergence', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin, Nbatch=nbFiles,
                         Np=nP, xdiv=ZDiv, quadArea=ZAqd )

    print(' *** Saving: '+cd_in+'/def_TOT_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_TOT_'+cfroot+'.npz', name='shear', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin, Nbatch=nbFiles,
                         Np=nP, xtot=Ztot, quadArea=ZAqd )


    
