#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir, environ
from glob import glob
import numpy as np
#from re import split
import mojito   as mjt

idebug=1
iplot=1

dir_in = '<HOME>/Nextcloud/data/mojitoNEW'

l_cst_bins = False ; rfexp_bin = 0.2

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

#do_scales = [ 10, 20, 40, 320, 640 ]
do_scales = [ 320, 640 ]

cfield = 'total'; cfld = 'tot'; cFLD = 'TOT'


if __name__ == '__main__':

    #if not len(argv) in [2]:
    #    print('Usage: '+argv[0]+' <file_divergence.npz>')
    #    exit(0)
    #cf_div_in = argv[1]
    #cf_shr_in = str.replace( cf_div_in, 'DIV', 'SHR' )


    dir_in = str.replace( dir_in, '<HOME>', environ.get('HOME') )

    print('\n *** Will find deformation files into: '+dir_in)


    Nscl = len(do_scales)


    xMean = np.zeros((Nscl,3))
    xVarc = np.zeros((Nscl,3))
    xSkew = np.zeros((Nscl,3))


    # Populating files:
    iscl = 0
    for res in do_scales:
        print('\n *** Resolution = '+str(res)+' km:')
        
        cc  = dir_in+'/rgps/def_'+cFLD+'_RGPS_dt*_'+str(res)+'km_????????-????????.npz'
        lst = np.sort( glob(cc) )
        if len(lst)!=1:
            print('ERROR: we do not have a single file!!! =>',cc)
        cf = lst[0]
    
        print(cf)

        with np.load(cf) as data:
            dtbin = data['dtbin']
            reskm = data['reskm_nmnl']
            corig = str(data['origin'])
            cperd = str(data['period'])        
            nbF   = int(data['Nbatch'])
            Ztot   =     data['xtot']

        (Ns,) = np.shape(Ztot) ; # sample size N
        #print('shape(Ztot) =',Ns)

        zm2   = Ztot*Ztot
        
        xMean[iscl,0] = np.mean(Ztot)        
        xVarc[iscl,0] = np.mean(zm2)
        xSkew[iscl,0] = np.mean(zm2*Ztot)
        
        print(' * Mean, Variance, Skewness =',xMean[iscl,0], xVarc[iscl,0], xSkew[iscl,0])


        iscl += 1
    exit(0)


    #cfroot = 'SCALING_'+corigin+'_dt'+cdtbin+'_'+str(reskm)+'km_'+cperiod   
    
