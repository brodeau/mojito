#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir, environ
from glob import glob
import numpy as np
from re import split
import mojito   as mjt

idebug=1
iplot=1

#dir_in = '<HOME>/Nextcloud/data/mojitoNEW'

l_cst_bins = False ; rfexp_bin = 0.2

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

cfield = 'total'; cfld = 'tot'; cFLD = 'TOT'
#cfield = 'shear'; cfld = 'shr'; cFLD = 'SHR'
#cfield = 'divergence'; cfld = 'div'; cFLD = 'DIV'

vORIGS = ['RGPS','BBM','EVP']


if __name__ == '__main__':

    if not len(argv)==3:
        print('Usage: '+argv[0]+' <dir_npz_in> <list_scales>')
        exit(0)
    dir_npz_in = argv[1]
    lst_scales = argv[2]
    do_scales = np.array( split(',',lst_scales), dtype=int )
    
    #dir_npz_in = str.replace( dir_npz_in, '<HOME>', environ.get('HOME') )

    print('\n *** Will find deformation files into: '+dir_npz_in)

    print('\n *** Scales to work with:', do_scales)

    Nscl = len(do_scales)

    xMQ = np.zeros((Nscl,3,3)) -9999. ; # Moments: [scale,origin,order]
    xAq = np.zeros((Nscl,3)  ) -9999. ; # Mean quad area: [scale,origin]
    reskm_actual = np.zeros((Nscl,3)  ) -9999. ; # actual resolution we're dealing with
    
    # Populating files:
    iscl = 0
    for res in do_scales:
        print('\n *** Resolution = '+str(res)+' km:')

        io = 0
        for corig in vORIGS:
            csdir = corig.lower()
            dirin = dir_npz_in+'/'+csdir
            if not path.exists(dirin):
                print('ERROR: directory "'+dirin+'" does not exist!'); exit(0)            
            cc  = dirin+'/def_'+cFLD+'_*'+corig+'*_dt*_'+str(res)+'km_????????-????????.npz'
            lst = np.sort( glob(cc) )
            if len(lst)!=1:
                print('ERROR: we do not have a single file!!! =>',cc); exit(0)
            cf = lst[0]

            print(cf)

            with np.load(cf) as data:
                dtbin = data['dtbin']
                reskm = data['reskm_nmnl']
                corig = str(data['origin'])
                cperd = str(data['period'])
                nbF   = int(data['Nbatch'])
                Ztot  =     data['x'+cfld]
                ZA    =     data['quadArea']

            (Ns,) = np.shape(Ztot) ; # sample size N
            #print('shape(Ztot) =',Ns)

            if cfield=='divergence':
                Ztot = np.abs(Ztot)

            zm2   = Ztot*Ztot

            xMQ[iscl,io,0], xMQ[iscl,io,1], xMQ[iscl,io,2] = np.mean(Ztot), np.mean(zm2), np.mean(zm2*Ztot)

            xAq[iscl,io] = np.mean(ZA)
            
            print(' * Mean, Variance, Skewness, mean quad area =',xMQ[iscl,io,0], xMQ[iscl,io,1], xMQ[iscl,io,2], xAq[iscl,io],'km^2')
            
            
            #if res==320 and io==0:
            #    for rr in Ztot:
            #        print(rr)
            #    print('')
            #    exit(0)


            del Ztot, ZA

            io += 1

        iscl += 1



    reskm_actual[:,:] = np.sqrt(xAq[:,:])
    
    ### Well I guess time for plot:


    

    if not path.exists('./figs'):
        mkdir('./figs')

    # Separate: 
    #cfroot = './figs/SCALING_'+cfield+'_mean_'+corig+'_dt'+str(dtbin)
    #kk = mjt.plotScalingDef( reskm_actual, xMQ[:,:,0], vORIGS, what='Mean', cfig=cfroot+'.svg' )
    #cfroot = './figs/SCALING_'+cfield+'_variance_'+corig+'_dt'+str(dtbin)
    #kk = mjt.plotScalingDef( reskm_actual, xMQ[:,:,1], vORIGS, what='Variance', cfig=cfroot+'.svg' )
    #cfroot = './figs/SCALING_'+cfield+'_skewness_'+corig+'_dt'+str(dtbin)
    #kk = mjt.plotScalingDef( reskm_actual, xMQ[:,:,2], vORIGS, what='Skewness', cfig=cfroot+'.svg' )

    cfroot = './figs/SCALING_'+cfield+'_'+corig+'_dt'+str(dtbin)
    kk = mjt.plot3ScalingDef( reskm_actual, xMQ, vORIGS, cfig=cfroot+'.svg' )
