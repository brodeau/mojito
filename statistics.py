#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
from glob import glob
import numpy as np
from re import split
#
from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=1
iplot=1


l_cst_bins = True
#l_cst_bins = False ; rfexp_bin = 0.1

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

# Conversion from s-1 to day-1:
rconv = 24.*3600.

# Range to consider:
min_div  = 0.001 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
max_div  = 5.     ; # day^-1
#
min_shr = 0.001 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
max_shr = 5. ; # day^-1


# About width of bins:
if l_cst_bins:
    wbin_div = 0.0005 ; # day^-1
    wbin_shr = 0.0005 ; # day^-1
else:
    wbin_min = 0.0005 ; # Narrowest bin width (for the smalles values of deformation)



### 
def constructCstBins( rmin, rmax, wdthB, name='divergence', iverbose=0 ):
    #
    nB = (rmax - rmin) / wdthB
    if not nB%1.==0.:
        print('ERROR [constructBins()]: "'+name+'" => nB is not an integer! nB =',nB); exit(0)
    nB = int(nB)
    #
    zbin_bounds = [  rmin + float(i)*wdthB for i in range(nB+1) ]
    zbin_bounds = np.round( zbin_bounds, 6 )
    zbin_center = [ 1.5*wdthB + float(i)*wdthB for i in range(nB) ]
    zbin_center = np.round( zbin_center, 6 )
    #
    if iverbose>0:
        print('\n * constructBins()]: we have '+str(nB)+' bins for the '+name+' !')
        print('     => zbin_bounds =',zbin_bounds,'\n')
        print('     => zbin_center =',zbin_center)

    # len(zbin_bounds) = nB + 1 !!!
    return nB, zbin_bounds, zbin_center


def constructExpBins( rfexp, rmin, rmax, wbmin, name='divergence', iverbose=0 ):
    #
    from math import exp
    #
    kc  = 0
    zwbin = []
    zacc  = rmin
    while not zacc>=rmax:
        rc = float(kc)
        wb = wbmin*exp(rfexp*rc)
        zwbin.append(wb)
        zacc = zacc+wb
        kc = kc+1

    zwbin = np.array(zwbin)
    nB    = len(zwbin)
    zbin_bounds = np.zeros(nB+1)
    zbin_center = np.zeros(nB)

    zbin_bounds[0] = rmin
    for jb in range(nB):
        zbin_bounds[jb+1] = zbin_bounds[jb] + zwbin[jb]
        zbin_center[jb]   = zbin_bounds[jb] + 0.5*zwbin[jb]

    if iverbose>0:
        for jb in range(nB):
            print( ' * [constructExpBins] Bin #',jb,' =>',zwbin[jb], zbin_bounds[jb],zbin_bounds[jb+1], zbin_center[jb] )

    print(' * [constructExpBins]: narrowest and widest bins:',zwbin[0],zwbin[-1])
    
    return nB, zbin_bounds, zbin_center



def computePDF( pBb, pBc, pX, cwhat='unknown', iverbose=0 ):
    '''
        * pBb, pBc: vectors for boundaries and center of bins
    '''
    (nBins,) = np.shape(pBc)
    if np.shape(pBb)!=(nBins+1,):
        print('ERROR [computePDF()]: `shape(pBb)!=(nBins+1,)` !!!')
        exit(0)
    
    (nbP,) = np.shape(pX)

    zxmin, zxmax = pBb[0], pBb[-1]
    zxrng = zxmax - zxmin
    zbwdth  = np.array( [ (pBb[i+1]-pBb[i]) for i in range(nBins) ] ) ; # width of bins...
    zbwdthN = zbwdth / (zxrng/nBins) ; # normalized width of bins... => should be only `1` when constant bins...
    #print( 'LOLO zbwdthN =',zbwdthN[:]) ;    exit(0)

    zPDF   = np.zeros(nBins)

    nPok = 0
    for iP in range(nbP):
    
        zX = pX[iP] ; # Yes! Absolute value of divergence, we want 1 single tail for the PDF...
    
        if zX>zxmax or zX<zxmin:
            if iverbose>0:
                print(' * WARNING [computePDF()]: excluding tiny or extreme value of '+cwhat+': ',zX,'day^-1')
            #
        else:
            jf  = np.argmin( np.abs( pBc - zX ) )
            lOk = ( zX>pBb[jf] and zX<=pBb[jf+1] )
            if not lOk:
                if zX<=pBb[jf]:   jf = jf-1
                if zX >pBb[jf+1]: jf = jf+1                
            if not ( zX>pBb[jf] and zX<=pBb[jf+1] ):
                print(' Binning error on divergence!')
                print('  => divergence =',zX)
                print('  => bounds =',pBb[jf],pBb[jf+1])
                exit(0)
            zPDF[jf] = zPDF[jf] + 1.
            
            nPok = nPok + 1; # another valid point

    # Normalization:
    
    #zPDF[:] =  zPDF[:]/float(nPok) ; # case constant bin width...
    
    #print(zbwdthN[0])
    #print('LOLO: BEFORE: integral of PDF =',np.sum(zPDF[:]*zbwdthN[:]))    
    #zPDF[:] = zPDF[:]*zbwdthN[:] / np.sum(zPDF[:]*zbwdthN[:])
    zPDF[:] = zPDF[:]*zbwdth[:]/zbwdthN[:] / np.sum(zPDF[:]*zbwdth[:]/zbwdthN[:])

    print('LOLO: AFTER: integral of PDF =',np.sum(zPDF[:]*zbwdthN[:]))
    #exit(0)
    
    return nPok, zPDF




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
        if kf==0: corigin = vf[1] ; # should be 'RGPS' or 'SI3' !
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
    cperiod = cdt1+'-'+cdt2

    
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

    cxtra = ''
    if l_cst_bins:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = constructCstBins( min_div, max_div, wbin_div, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = constructCstBins( min_shr, max_shr, wbin_shr, name='shear', iverbose=idebug )
        
    else:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = constructExpBins( rfexp_bin, min_div, max_div, wbin_min, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = constructExpBins( rfexp_bin, min_shr, max_shr, wbin_min, name='shear', iverbose=idebug )
        cxtra = '_incB'
    

    nPd, PDF_div = computePDF( xbin_bounds_div, xbin_center_div, np.abs(Zdiv), cwhat='divergence' )

    nPs, PDF_shr = computePDF( xbin_bounds_shr, xbin_center_shr,        Zshr , cwhat='shear' )

        

    
    #print(' Integral value of PDF_div =',np.sum(PDF_div[:]*zbwdth[:] ) )
    #exit(0)
        
    #else:
          #

    #zd = np.sum(PDF_div[:]*zbwdth[:]
        
    #PDF_div[:] = PDF_div[:]*zbwdth[:]/np.sum(zbwdth)
    #zbwdth = np.array( [ xbin_bounds_shr[i+1]-xbin_bounds_shr[i] for i in range(nBinsS) ] )        
    #PDF_shr[:] = PDF_shr[:]*zbwdth[:]/np.sum(zbwdth)

    #PDF_div[:] = PDF_div[:] / np.sum(PDF_div[:]*zbwdth[:])
    


    cfroot = 'PDF_'+corigin+'_'+cperiod



    
    # Saving in `npz` files:
    np.savez_compressed( './npz/'+cfroot+'_divergence.npz', name='divergence', origin=corigin, period=cperiod, wbin=0,
                         Np=nPd, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_div )
    
    np.savez_compressed( './npz/'+cfroot+'_shear.npz',      name='shear',      origin=corigin, period=cperiod, wbin=0,
                         Np=nPs, xbin_bounds=xbin_bounds_shr, xbin_center=xbin_center_shr, PDF=PDF_shr )


    if iplot>0:
        cdir = './figs'
        if not path.exists(cdir): mkdir(cdir)
    
            
        kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, Np=nPd, name='Divergence', cfig=cdir+'/loglog'+cfroot+'_divergence'+cxtra+'.png',
                            wbin=0, title=corigin, period=cperiod )
        
        kk = mjt.LogPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPd, name='Shear', cfig=cdir+'/loglog'+cfroot+'_shear'+cxtra+'.png',
                            wbin=0, title=corigin, period=cperiod )
    
        xdiv_rng=[0.001,0.1] ; # x-range we want on the x-axis of the plot
        xshr_rng=[0.001,0.1] ; # x-range we want on the x-axis of the plot
        
        kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, Np=nPd, name='Divergence', cfig=cdir+'/'+cfroot+'_divergence'+cxtra+'.png',
                             xrng=xdiv_rng, wbin=0, title=corigin, period=cperiod )
        
        kk = mjt.PlotPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPs, name='Shear', cfig=cdir+'/'+cfroot+'_shear'+cxtra+'.png',
                             xrng=xshr_rng, wbin=0, title=corigin, period=cperiod )
    
