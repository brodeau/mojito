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

    if not len(argv) in [3,5]:
        print('Usage: '+argv[0]+' <directory_input_npz_files> <dtbin_h> <creskm> <string_id_origin>')
        print('   or: '+argv[0]+' <directory_input_npz_files> <file_prefix>')
        exit(0)
    
    lPrefix = (len(argv)==3)
    cd_in  = argv[1]

    if lPrefix:
        cprfx = argv[2]
        listnpz = np.sort( glob(cd_in+'/'+cprfx+'*.npz') )
        cdtbin, creskm = '', ''
        dtbin = 0
    else:
        cdtbin = argv[2]
        creskm = argv[3]
        cidorg = argv[4]
        listnpz = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'*_dt'+cdtbin+'*'+creskm+'km.npz') )

    # Polpulating deformation files available:    
    nbFiles = len(listnpz)
    print('\n *** We found '+str(nbFiles)+' deformation files into '+cd_in+' !')

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

    if not lPrefix:
        dtbin=int(cdtbin)
        if str(reskm) != creskm:
            print('ERROR: spatial scale (km) passed as argument does not match that found in deformation files!')
            exit(0)
    
    # Now that we know the total number of points we can allocate and fill arrays for divergence and shear    
    Zshr = np.zeros(nP)
    ZDiv = np.zeros(nP)
    Ztot = np.zeros(nP)
    
    jP = 0 ; # Counter from 0 to nP-1
    kf = 0 ; # Counter for files, 0 to nbFiles-1
    for ff in listnpz:
        jPe = jP + kNbPoints[kf]
        with np.load(ff) as data:
            zdiv  =      data['divergence']
            zshr  =      data['shear']
            ztot  =      data['total']
        #
        ZDiv[jP:jPe] = rconv*zdiv ; # day^-1
        Zshr[jP:jPe] = rconv*zshr ; # day^-1
        Ztot[jP:jPe] = rconv*ztot ; # day^-1
        #
        jP = jPe
        kf = kf+1
    

    if lPrefix:
        cfroot = corigin+'_'+str(reskm)+'km_'+cperiod
    else:
        cfroot = corigin+'_dt'+cdtbin+'_'+str(reskm)+'km_'+cperiod   

    
    # Save it for the scaling and PDFs to come ...
    np.savez_compressed( cd_in+'/def_SHR_'+cfroot+'.npz', name='shear', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin, Nbatch=nbFiles,
                         Np=nP, xshr=Zshr )

    np.savez_compressed( cd_in+'/def_DIV_'+cfroot+'.npz', name='divergence', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin, Nbatch=nbFiles,
                         Np=nP, xdiv=ZDiv )

    np.savez_compressed( cd_in+'/def_TOT_'+cfroot+'.npz', name='shear', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin, Nbatch=nbFiles,
                         Np=nP, xtot=Ztot )


    
