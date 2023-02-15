#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
from glob import glob
import numpy as np
from re import split
from math import ceil,floor

#from scipy.spatial import Delaunay
import matplotlib.pyplot as plt


from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=1

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

# Conversion from s-1 to day-1:
rconv = 24.*3600.


# Bin widths for pdfs
wbin_div = 0.0005 ; # day^-1
min_div  = 0.001 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
max_div  = 5.     ; # day^-1

wbin_shr = 0.0005 ; # day^-1
min_shr = 0.001 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
max_shr = 5. ; # day^-1




### 
def constructBins( rmin, rmax, wdthB, name='divergence', iverbose=0 ):
    #
    nB = (rmax - rmin) / wdthB
    if not nB%1.==0.:
        print('ERROR [constructBins()]: "'+name+'" => nB is not an integer! nB =',nB); exit(0)
    nB = int(nB)
    #
    zbin_bounds = [  float(i+1)*wdthB for i in range(nB+1) ]
    zbin_bounds = np.round( zbin_bounds, 6 )
    zbin_center = [ 1.5*wdthB + float(i)*wdthB for i in range(nB) ]
    zbin_center = np.round( zbin_center, 6 )
    #
    if iverbose>0:
        print('\n * constructBins()]: we have '+str(nB)+' bins for the '+name+' !')
        print('     => zbin_bounds =',zbin_bounds,'\n')
        print('     => zbin_center =',zbin_center)

    return nB, zbin_bounds, zbin_center






if __name__ == '__main__':

    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <directory_input_npz_files>')
        exit(0)
    cd_in = argv[1]


    # Polpulating deformation files available:
    listnpz = np.sort( glob(cd_in+'/'+cprefixIn+'*.npz') )
    nbFiles = len(listnpz)
    print('\n *** We found '+str(nbFiles)+' deformation files into '+cd_in+' !')


    kStreamName = np.zeros(nbFiles, dtype='U4')
    kiDate      = np.zeros(nbFiles, dtype=int ) ; # date in epoch time at which deformations were calculated
    kNbPoints   = np.zeros(nbFiles, dtype=int ) ; # number of points in file

    list_date = []
    kf = 0
    for ff in listnpz:
        print('\n  # File: '+ff)
        fb = path.basename(ff)
        vf = split('_',fb)
        print(vf)
        #
        if kf==0: cname = vf[1] ; # should be 'RGPS' or 'SI3' !
        list_date.append(split('-',vf[2])[0])
        #
        kStreamName[kf] = vf[1]
        #
        with np.load(ff) as data:
            rdate = int( data['time'] )
            nPnts =      data['Npoints']
    
        kiDate[kf] = rdate
        kNbPoints[kf] = nPnts
    
        print('   * Stream: '+kStreamName[kf] )
        print('   * Date = ',epoch2clock(kiDate[kf]))
        print('   * Nb. of points = ',kNbPoints[kf] )
            
        kf = kf+1
    
    print('\n')
    
    #print('  ==> list of streams:', kStreamName[:])
    
    nP = np.sum(kNbPoints)
    print('  ==> Total number of points:', nP)


    print('  ==> list of dates:', list_date[:])

    cdt1, cdt2 = list_date[0],list_date[-1]
    
    # Now that we know the total number of points we can allocate and fill arrays for divergence and shear
    Zdiv = np.zeros(nP)
    Zshr = np.zeros(nP)
    
    jP = 0 ; # Counter from 0 to nP-1
    kf = 0 ; # Counter for files, 0 to nbFiles-1
    for ff in listnpz:
        jPe = jP + kNbPoints[kf]
        with np.load(ff) as data:
            zdiv  =      data['divergence']
            zshr  =      data['shear']
        #
        Zdiv[jP:jPe] = rconv*zdiv ; # day^-1
        Zshr[jP:jPe] = rconv*zshr ; # day^-1
        #
        jP = jPe
        kf = kf+1


    # For the divergence
    nBinsD, xbin_bounds_div, xbin_center_div = constructBins( min_div, max_div, wbin_div, name='divergence', iverbose=idebug )
        
    # For the shear:
    nBinsS, xbin_bounds_shr, xbin_center_shr = constructBins( min_shr, max_shr, wbin_shr, name='shear', iverbose=idebug )
    
    
    PDF_div = np.zeros(nBinsD)
    PDF_shr = np.zeros(nBinsS)
    
    nPd = nP
    nPs = nP
    
    for iP in range(nP):
    
        rdiv = abs(Zdiv[iP]) ; # Yes! Absolute value of divergence, we want 1 single tail for the PDF...
    
        if rdiv > max_div:
            nPd = nPd - 1
            print(' * WARNING: excluding extreme value of Divergence: ',rdiv,'day^-1')
        elif rdiv <= min_div:
            nPd = nPd - 1
            #
        else:
            jf = np.argmin( np.abs( xbin_center_div - rdiv ) )    
            if not ( rdiv>xbin_bounds_div[jf] and rdiv<=xbin_bounds_div[jf+1] ):
                print(' Binning error on divergence!')
                print('  => divergence =',rdiv)
                print('  => bounds =',xbin_bounds_div[jf],xbin_bounds_div[jf+1])
                exit(0)
            PDF_div[jf] = PDF_div[jf]+1
    
        rshr = Zshr[iP]
        if rshr> max_shr or rshr<0.:
            nPs = nPs - 1
            print(' * WARNING: excluding extreme value of Shear: ',rshr,'day^-1')        
        elif rshr <= min_shr:
            nPd = nPd - 1
        else:
            jf = np.argmin( np.abs( xbin_center_shr - rshr ) )    
            if not ( rshr>=xbin_bounds_shr[jf] and rshr<xbin_bounds_shr[jf+1] ):
                print(' Binning error on shear!')
                print('  => shear =',rshr)
                print('  => bounds =',xbin_bounds_shr[jf],xbin_bounds_shr[jf+1])
                exit(0)
            PDF_shr[jf] = PDF_shr[jf]+1
    
    
    
    PDF_div[:] = PDF_div[:]/float(nPd)
    PDF_shr[:] = PDF_shr[:]/float(nPs)



    cfroot = 'PDF_'+cname+'_'+cdt1+'-'+cdt2



    
    # Saving in `npz` files:
    np.savez_compressed( './npz/'+cfroot+'_divergence.npz', name='divergence', Np=nPd, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_div )
    np.savez_compressed( './npz/'+cfroot+'_shear.npz',      name='shear',      Np=nPs, xbin_bounds=xbin_bounds_shr, xbin_center=xbin_center_shr, PDF=PDF_shr )

        
    kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, Np=nPd, name='Divergence', cfig='loglog'+cfroot+'_divergence.svg',
                        wbin=wbin_div, title=cname, period=cdt1+' - '+cdt2 )
    
    kk = mjt.LogPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPd, name='Shear', cfig='loglog'+cfroot+'_shear.svg',
                        wbin=wbin_shr, title=cname, period=cdt1+' - '+cdt2 )


    xdiv_rng=[min_div,0.1] ; # x-range we want on the x-axis of the plot
    xshr_rng=[min_shr,0.1] ; # x-range we want on the x-axis of the plot
    
    kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, Np=nPd, name='Divergence', cfig=cfroot+'_divergence.svg',
                         xrng=xdiv_rng, wbin=wbin_div, title=cname, period=cdt1+' - '+cdt2 )
    
    kk = mjt.PlotPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPs, name='Shear', cfig=cfroot+'_shear.svg',
                         xrng=xshr_rng, wbin=wbin_shr, title=cname, period=cdt1+' - '+cdt2 )
    
