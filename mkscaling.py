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

do_scales = np.array([ 10, 20, 40, 80, 160, 320, 640 ], dtype=int)

cfield = 'total'; cfld = 'tot'; cFLD = 'TOT'
#cfield = 'shear'; cfld = 'shr'; cFLD = 'SHR'
#cfield = 'divergence'; cfld = 'div'; cFLD = 'DIV'

vORIGS = ['RGPS','BBM','EVP']


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

        io = 0
        for corig in vORIGS:
            csdir = corig.lower()

            cc  = dir_in+'/'+csdir+'/def_'+cFLD+'_*'+corig+'*_dt*_'+str(res)+'km_????????-????????.npz'
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
                Ztot   =     data['x'+cfld]

            (Ns,) = np.shape(Ztot) ; # sample size N
            #print('shape(Ztot) =',Ns)

            if cfield=='divergence':
                Ztot = np.abs(Ztot)

            zm2   = Ztot*Ztot

            xMean[iscl,io] = np.mean(Ztot)
            xVarc[iscl,io] = np.mean(zm2)
            xSkew[iscl,io] = np.mean(zm2*Ztot)

            print(' * Mean, Variance, Skewness =',xMean[iscl,io], xVarc[iscl,io], xSkew[iscl,io])


            #if res==320 and io==0:
            #    for rr in Ztot:
            #        print(rr)
            #    print('')
            #    exit(0)


            del Ztot

            io += 1

        iscl += 1

    ### Well I guess time for plot:

    cfroot = './figs/SCALING_'+cfield+'_mean_'+corig+'_dt'+str(dtbin)
    kk = mjt.plotScalingDef( do_scales, xMean, vORIGS, what='Mean', cfig=cfroot+'.png' )

    cfroot = './figs/SCALING_'+cfield+'_variance_'+corig+'_dt'+str(dtbin)
    kk = mjt.plotScalingDef( do_scales, xVarc, vORIGS, what='Variance', cfig=cfroot+'.png' )

    cfroot = './figs/SCALING_'+cfield+'_skewness_'+corig+'_dt'+str(dtbin)
    kk = mjt.plotScalingDef( do_scales, xSkew, vORIGS, what='Skewness', cfig=cfroot+'.png' )
