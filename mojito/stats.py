#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

#from sys import argv, exit
#from os import path, mkdir
#from glob import glob
import numpy as np
#from re import split
#import mojito   as mjt
from .util import epoch2clock as e2c

idebug=1
#iplot=1


def constructCstBins( rmin, rmax, wdthB, name='unknown', iverbose=0 ):
    #
    nB = (rmax - rmin) / wdthB
    if not nB%1.==0.:
        print('ERROR [constructBins()]: "'+name+'" => nB is not an integer! nB =',nB); exit(0)
    nB = int(nB)
    #
    zbin_bounds = [ rmin + float(i)*wdthB                 for i in range(nB+1) ]
    zbin_center = [ 0.5*(zbin_bounds[i]+zbin_bounds[i+1]) for i in range(nB)   ]    
    zbin_bounds = np.round( zbin_bounds, 6 )
    zbin_center = np.round( zbin_center, 6 )
    #
    if iverbose>0:
        for jb in range(nB):
            print( ' * [constructBins] Bin #',jb,' =>', zbin_bounds[jb],zbin_bounds[jb+1], zbin_center[jb] )
        
    # len(zbin_bounds) = nB + 1 !!!
    return nB, zbin_bounds, zbin_center


def constructExpBins( rfexp, rmin, rmax, wbmin, name='unknown', iverbose=0 ):
    #
    from math import exp
    #
    kc  = 0
    zwbin = []
    zacc  = rmin
    while not zacc>=rmax:
        rc = float(kc)
        wb = wbmin*exp(rfexp*rc)
        #wb = wbmin*10.**(rfexp*rc)
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



def computePDF( pBb, pBc, pX, cwhat='unknown', return_cleaned=False, iverbose=0 ):
    '''
        * pBb, pBc: vectors for boundaries and center of bins
    '''
    (nBins,) = np.shape(pBc)
    if np.shape(pBb)!=(nBins+1,):
        print('ERROR [computePDF()]: `shape(pBb)!=(nBins+1,)` !!!')
        exit(0)
    
    (nbP,) = np.shape(pX)
    #
    zxmin, zxmax = pBb[0], pBb[-1]
    zxrng = zxmax - zxmin
    zbW  = np.array( [ (pBb[i+1]-pBb[i]) for i in range(nBins) ] ) ; # width of bins...
    
    lvbw = (  np.sum(np.abs(zbW[1:]-zbW[:-1])) > 1.e-11 ) ; # Constant or variable-width bins?
    if lvbw:
        if iverbose>0: print(' * [computePDF()]: variable-width bins for '+cwhat+'!')
        #zscal   = zbW[0]/zbW[:]
        zscal   = 1./zbW[:]
    #
    if return_cleaned:
        zXclean = []
    
    zPDF = np.zeros(nBins)
    nPok = 0
    for iP in range(nbP):
    
        zX = pX[iP]
    
        if zX>zxmax or zX<zxmin:
            if iverbose>1:
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
            #
            if return_cleaned:
                zXclean.append(zX)

            

    # Normalization:
    if lvbw:
        zPDF[:] = zPDF[:]*zscal[:]
    #
    zPDF[:] = zPDF[:]/float(nPok)

    if return_cleaned:
        return nPok, zPDF, np.array(zXclean)
    else:
        return nPok, zPDF





def Construct90P( ifile, vdates_batch, pdates, pdef, pp=90, Nmin=1000, iverbose=0  ):
    '''
        * ifile: number of the file treated
        *  pdef must be in days^-1
        
    '''
    lmaskArrays = True
    (nbF,) = np.shape(vdates_batch)
    
    Z90P = np.zeros( nbF )
    zdat = np.zeros( nbF, dtype=int )
    imsk = np.zeros( nbF, dtype='i1' )

    zdates, zdef = pdates.copy(), pdef.copy()

    # Need to remove erroneous extremely small values:
    from mojito import config as cfg
    kk = cfg.updateConfig4Scale( 10, mode='rgps', ltalk=False )
    print(' * [Construct90P()] => uses `rc_MinDef` =',cfg.rc_MinDef,'to clean all data / rgps, regardless of origin!!!')        
    (idxKeep,) = np.where(pdef>cfg.rc_MinDef)
    zdef = pdef[idxKeep]
    zdates = pdates[idxKeep]

    
    ic = 0
    for jd in vdates_batch:
        if iverbose>0: print('\n *** file #'+str(ifile)+''+e2c(jd))

        (idxDate,) = np.where( zdates == jd )
        nV = len(idxDate)
        if iverbose>0: print('         => '+str(nV)+' deformation for this date....')
        
        if nV>=Nmin:
            # sample size must be large enough
            ztmp = zdef[idxDate]
            Z90P[ic] = np.percentile(ztmp, pp)
            zdat[ic] =  jd
            imsk[ic] = 1
            
        ic+=1

    if lmaskArrays:
        Z90Pm = np.ma.masked_where(imsk!=1, Z90P)
        zdatm = np.ma.masked_where(imsk!=1, zdat)
    
        zt = np.ma.MaskedArray.compressed( zdatm )
        zx = zt[1:] - zt[:-1]
        if np.any(zx<=0):
            print('ERROR [Construct90P]: `zdatm` is not increasing!'); exit(0)

        return zt, np.ma.MaskedArray.compressed( Z90Pm )
    else:
        (idxM,) = np.where(imsk==0)
        Z90P[idxM] = 0.

        return zdat, Z90P





def Construct90P_old( ifile, vdates_batch, pdates, pdef, Nmin=1000 ):
    '''
        * ifile: number of the file treated
        
    '''
    lmaskArrays = True
    (nbF,) = np.shape(vdates_batch)
    
    Z90P = np.zeros( nbF )
    zdat = np.zeros( nbF, dtype=int )
    imsk = np.zeros( nbF, dtype='i1' )

    zdates, zdef = pdates.copy(), pdef.copy()

    ic = 0
    for jd in vdates_batch:
        print('\n *** file #'+str(ifile)+''+e2c(jd))

        (idxDate,) = np.where( zdates == jd )
        nV = len(idxDate)
        print('         => '+str(nV)+' deformation for this date....')
        
        ztmp = zdef[idxDate]
        zdfw = np.sort(ztmp) ; # Sorted in increasing order!
        #for rr in zdfw:
        #    print(rr)
        #exit(0)
        
        if zdfw.shape != (nV,):
            print('ERROR [Construct90P_old()]: problem #1'); exit(0)

        ri90 = 0.9*float(nV)
        i90f = int( floor( ri90 ) ) - 1   ; # `-1` because C-indexing...
        rw = ri90 - float(i90f)           ; # weight for interpolation

        if nV>=Nmin and i90f<nV-2:
            # sample size must be large enough
            # and: otherwize we are too close to the end of the series, this probably a bad batch!
            z90 = (1.-rw)*zdfw[i90f] + rw*zdfw[i90f+1]
            
            Z90P[ic] = z90
            zdat[ic] =  jd
            imsk[ic] = 1
            if z90>=Nmin:
                imsk[ic] = 1            
            #print('         => z90 =', z90, '(mask =',imsk[ic],')')
            
        ic+=1

    if lmaskArrays:
        Z90Pm = np.ma.masked_where(imsk!=1, Z90P)
        zdatm = np.ma.masked_where(imsk!=1, zdat)
    
        zt = np.ma.MaskedArray.compressed( zdatm )
        zx = zt[1:] - zt[:-1]
        if np.any(zx<=0):
            print('ERROR [Construct90P_old]: `zdatm` is not increasing!'); exit(0)
        print('')        
        return zt, np.ma.MaskedArray.compressed( Z90Pm )
    else:
        (idxM,) = np.where(imsk==0)
        Z90P[idxM] = 0.

        return zdat, Z90P

    

