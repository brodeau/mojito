#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir, environ
from glob import glob
import numpy as np
from re import split
import mojito   as mjt
from mojito import config as cfg

idebug=1
iplot=1

lAdaptMinDef = False

#dir_in = '<HOME>/Nextcloud/data/mojitoNEW'

l_cst_bins = False ; rfexp_bin = 0.2

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

cfield = 'total'; cfld = 'tot'; cFLD = 'TOT'
#cfield = 'shear'; cfld = 'shr'; cFLD = 'SHR'
#cfield = 'divergence'; cfld = 'div'; cFLD = 'DIV'

lOnlyRGPS = False

if __name__ == '__main__':

    k1 = cfg.initialize( mode='model' )
    
    if not len(argv)==4:
        print('Usage: '+argv[0]+' <dir_npz_in> <list_exp> <list_scales>')
        exit(0)
    dir_npz_in = argv[1]
    lst_exp    = argv[2]
    lst_scales = argv[3]
    
    do_scales = np.array( split(',',lst_scales), dtype=int  )
    vORIGS    = np.array( split(',',lst_exp),    dtype='U16')



    #dir_npz_in = str.replace( dir_npz_in, '<HOME>', environ.get('HOME') )

    print('\n *** Will find deformation files into: '+dir_npz_in)

    print('\n *** Scales to work with:', do_scales)

    Nscl = len(do_scales)
    Norg = len(vORIGS)

    cf  = np.zeros((Nscl,Norg), dtype='U256' )
    Nbp = np.zeros((Nscl,Norg), dtype=int    )
    xMQ = np.zeros((Nscl,Norg,3)) -9999. ; # Moments: [scale,origin,order]
    xAq = np.zeros((Nscl,Norg)  ) -9999. ; # Mean quad area: [scale,origin]
    reskm_actual = np.zeros((Nscl,Norg)  ) -9999. ; # actual resolution we're dealing with



    zmin_tot_def = np.zeros(Nscl) ; # minimum deformation allowed for a givn scale


    k2 = cfg.updateConfig4Scale( 10, mode='model' )
    zmin_tot_def[:] = cfg.rc_tot_min

    if lAdaptMinDef:
        iscl=0
        for res in do_scales:
            if res==20:
                zmin_tot_def[iscl] = 0.8e-4        
            if res==40:
                zmin_tot_def[iscl] = 1.5e-4        
            if res==80:
                zmin_tot_def[iscl] = 2.e-4        
            if res==160:
                zmin_tot_def[iscl] = 3.e-4        
            if res==320:
                zmin_tot_def[iscl] = 3.5e-4        
            if res==640:
                zmin_tot_def[iscl] = 1.5e-3
            iscl+=1

    print(zmin_tot_def[:])
    #exit(0)
        
    ### Save the name of files to read
    ### and get the number of points...
    iscl = 0
    for res in do_scales:
        print('\n *** Resolution = '+str(res)+' km:')

        io = 0
        for corig in vORIGS:
            csdir = corig.lower()
            dirin = dir_npz_in+'/'+csdir
            if not path.exists(dirin):
                print('ERROR: directory "'+dirin+'" does not exist!'); exit(0)            
            cc  = dirin+'/def_'+cFLD+'_*'+corig+'*_'+str(res)+'km_????????-????????.npz'
            lst = np.sort( glob(cc) )
            if len(lst)>1: print('ERROR: we have more than 1 file!!! =>',cc); exit(0)
            if len(lst)<1: print('ERROR: we do not have any file!!! =>',cc); exit(0)
            cf[iscl,io] = lst[0]
            cfi = cf[iscl,io]
            print(cfi)

            print('  => opening: '+cfi)
            with np.load(cfi) as data:
                Ztot  =     data['x'+cfld]

                
            (idxOK,) = np.where( np.abs(Ztot)>zmin_tot_def[iscl] )
            Nl = len(idxOK)
            Nbp[iscl,io] = Nl

            del Ztot, idxOK

            io += 1

        iscl += 1

    ###

    #print(' ***  cf =',cf)
    print(' *** Nbp =', Nbp)

    Nl_max = np.max(Nbp)
    print(' *** Biggest sample: ',str(Nl_max),'points!')

    xXQ = np.zeros((Nl_max,Nscl,Norg,3)) -9999. ; # Moments: [points,scale,origin,order]
    xXS = np.zeros((Nl_max,Nscl,Norg)) -9999. ; # Moments: [points,scale,origin]

    
    
    iscl = 0
    for res in do_scales:
        print('\n *** Resolution = '+str(res)+' km:')

        io = 0
        for corig in vORIGS:
            cfi = cf[iscl,io]
            
            print('  => opening: '+cfi)
            with np.load(cfi) as data:
                dtbin = data['dtbin']
                reskm = data['reskm_nmnl']
                corig = str(data['origin'])
                cperd = str(data['period'])
                nbF   = int(data['Nbatch'])
                Ztot  =     data['x'+cfld]
                ZA    =     data['quadArea']


            if cfield=='divergence':
                Ztot = np.abs(Ztot)

            # Removing bad tiny values:
            (idxOK,) = np.where( Ztot>zmin_tot_def[iscl] )
            Nl = len(idxOK) ; # sample size N
            Ztot = Ztot[idxOK]
            ZA   =   ZA[idxOK]
            
            zm2   = Ztot*Ztot
            zm3   = zm2*Ztot

            xMQ[iscl,io,0], xMQ[iscl,io,1], xMQ[iscl,io,2] = np.mean(Ztot), np.mean(zm2), np.mean(zm3)

            xAq[iscl,io] = np.mean(ZA)

            xXQ[:Nl,iscl,io,0], xXQ[:Nl,iscl,io,1], xXQ[:Nl,iscl,io,2] = Ztot[:], zm2[:], zm3[:]

            xXS[:Nl,iscl,io] = np.sqrt( ZA[:] )
            
            
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


    xXQ = np.ma.masked_where(xXQ<-9000., xXQ)


    Naxis = 0
    
    cfroot = './figs/SCALING_'+cfield+'_'+corig+'_dt'+str(dtbin)
    kk = mjt.plot3ScalingDef( reskm_actual, xMQ, vORIGS, pXQ=xXQ, pXS=xXS, cfig=cfroot+'.png',
                              lOnlyObs=lOnlyRGPS, lShowScat=True, Naxis=Naxis )

    # Separate: 
    cfroot = './figs/0scaling_'+cfield+'_q1_mean_'+corig+'_dt'+str(dtbin)
    kk = mjt.plotScalingDef( reskm_actual, xMQ[:,:,0], vORIGS, what='Mean', cfig=cfroot+'.png' )
    cfroot = './figs/0scaling_'+cfield+'_q2_variance_'+corig+'_dt'+str(dtbin)
    kk = mjt.plotScalingDef( reskm_actual, xMQ[:,:,1], vORIGS, what='Variance', cfig=cfroot+'.png' )
    cfroot = './figs/0scaling_'+cfield+'_q3_skewness_'+corig+'_dt'+str(dtbin)
    kk = mjt.plotScalingDef( reskm_actual, xMQ[:,:,2], vORIGS, what='Skewness', cfig=cfroot+'.png' )
    
