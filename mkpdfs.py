#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
#from glob import glob
import numpy as np
#from re import split
import mojito   as mjt

idebug=1
iplot=1

l_cst_bins = False ; rfexp_bin = 0.2

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

max_div = 1.5 ; # day^-1
max_shr = 1.5 ; # day^-1
max_tot = 1.5 ; # day^-1

# About width of bins:
if l_cst_bins:
    wbin_div = 0.001 ; # day^-1
    wbin_shr = 0.001 ; # day^-1
    wbin_tot = 0.001 ; # day^-1
else:
    wVbin_min = 0.0005 ; # Narrowest bin width (for the smalles values of deformation)

l_add_gaussian = False
    

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



if __name__ == '__main__':

    #if not len(argv) in [3,5]:
    #    print('Usage: '+argv[0]+' <directory_input_npz_files> <dtbin_h> <creskm> <string_id_origin>')
    #    print('   or: '+argv[0]+' <directory_input_npz_files> <file_prefix>')
    #    exit(0)
    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <file_divergence.npz>')
        exit(0)

    cf_div_in = argv[1]
    cf_shr_in = str.replace( cf_div_in, 'DIV', 'SHR' )
    cf_tot_in = str.replace( cf_div_in, 'DIV', 'TOT' )

    cd_in = path.dirname(cf_div_in)
    
    with np.load(cf_div_in) as data:
        dtbin1 = data['dtbin']
        reskm1 = data['reskm_nmnl']
        corig1 = str(data['origin'])        
        cperd1 = str(data['period'])
        nbF1 = int(data['Nbatch'])
        ZDiv = data['xdiv']

    with np.load(cf_shr_in) as data:
        dtbin2 = data['dtbin']
        reskm2 = data['reskm_nmnl']
        corig2 = str(data['origin'])
        cperd2 = str(data['period'])        
        nbF2   = int(data['Nbatch'])
        Zshr   =     data['xshr']

    with np.load(cf_tot_in) as data:
        dtbin3 = data['dtbin']
        reskm3 = data['reskm_nmnl']
        corig3 = str(data['origin'])
        cperd3 = str(data['period'])        
        nbF3   = int(data['Nbatch'])
        Ztot   =     data['xtot']

        
    # Control agreement:
    if dtbin1!=dtbin2 or dtbin1!=dtbin3:
        print('ERROR: `dtbin1!=dtbin2` !!!'); exit(0)
    if reskm1!=reskm2 or reskm1!=reskm3:
        print('ERROR: `reskm1!=reskm2` !!!'); exit(0)
    if corig1!=corig2 or corig1!=corig3:
        print('ERROR: `corig1!=corig2` !!!'); exit(0)
    if cperd1!=cperd2 or cperd1!=cperd3:
        print('ERROR: `cperd1!=cperd2` !!!'); exit(0)
    if nbF1!=nbF2 or nbF1!=nbF3:
        print('ERROR: `nbF1!=nbF2` !!!'); exit(0)

    dtbin  = dtbin1
    cdtbin = str(dtbin1)
    reskm  = reskm1
    creskm = str(reskm1)
    corigin = corig1
    cperiod = cperd1

    # For large scales, we must increase the size of bins:
    min_div, min_shr, min_tot = 0.003, 0.003, 0.003 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
    if   reskm >= 15. and reskm < 35.:
        min_div, min_shr, min_tot = 0.001, 0.001
    elif   reskm >= 35. and reskm < 50.:
        min_div, min_shr, min_tot = 5.e-4, 5.e-4, 5.e-4
        wVbin_min = 1.5*wVbin_min        
        rfexp_bin = rfexp_bin*1.2
    elif reskm >= 50. and reskm < 100.:
        min_div, min_shr, min_tot = 1.e-4, 1.e-4, 1.e-4
        wVbin_min = 2.*wVbin_min
        rfexp_bin = rfexp_bin*1.5
    elif reskm >= 100. and reskm < 200.:
        min_div, min_shr, min_tot = 1.e-5, 1.e-5, 1.e-5
        wVbin_min = 4.*wVbin_min
        rfexp_bin = rfexp_bin*2.
    elif reskm >= 200. and reskm < 400.:
        min_div, min_shr, min_tot = 1.e-5, 1.e-5, 1.e-5
        wVbin_min = 6.*wVbin_min
        rfexp_bin = rfexp_bin*2.5
    elif reskm >= 400. and reskm < 700.:
        min_div, min_shr, min_tot = 1.e-5, 1.e-5, 1.e-5
        wVbin_min = 8.*wVbin_min
        rfexp_bin = rfexp_bin*2.5
    if reskm >= 35.:
        print('\n * Increased the width of bins and exp growth!!! Because large scale! wVbin_min, rfexp_bin =',wVbin_min, rfexp_bin)
            
    print(ZDiv)

    
    cxtra = ''
    if l_cst_bins:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = constructCstBins( min_div, max_div, wbin_div, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = constructCstBins( min_shr, max_shr, wbin_shr, name='shear',      iverbose=idebug )        
        # For the total deformation:
        nBinsS, xbin_bounds_tot, xbin_center_tot = constructCstBins( min_tot, max_tot, wbin_tot, name='shear',      iverbose=idebug )
        
    else:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = constructExpBins( rfexp_bin, min_div, max_div, wVbin_min, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = constructExpBins( rfexp_bin, min_shr, max_shr, wVbin_min, name='shear',      iverbose=idebug )
        # For the total deformation:
        nBinsS, xbin_bounds_tot, xbin_center_tot = constructExpBins( rfexp_bin, min_tot, max_tot, wVbin_min, name='total',      iverbose=idebug )
        cxtra = '_incB'



    # Signed divergence:
    (idxN,) = np.where(ZDiv<-min_div)
    nPn  = len(idxN)
    Zcnv = np.zeros(nPn)
    Zcnv[:] = - ZDiv[idxN] ; # We want it to be positive for the graph...
    
    (idxP,) = np.where(ZDiv> min_div)
    nPp  = len(idxP)
    Zdiv = np.zeros(nPp)
    Zdiv[:] = ZDiv[idxP]


    if l_add_gaussian:
        nPs, PDF_shr, ZshrClean = computePDF( xbin_bounds_shr, xbin_center_shr, Zshr, cwhat='shear', return_cleaned=True   ,  iverbose=idebug )
    else:
        nPs, PDF_shr = computePDF( xbin_bounds_shr, xbin_center_shr,    Zshr,  cwhat='shear',      iverbose=idebug )

    nPs, PDF_tot = computePDF( xbin_bounds_tot, xbin_center_tot,        Ztot,  cwhat='total',       iverbose=idebug )    
    nPD, PDF_Div = computePDF( xbin_bounds_div, xbin_center_div, np.abs(ZDiv), cwhat='Divergence',  iverbose=idebug )
    nPd, PDF_div = computePDF( xbin_bounds_div, xbin_center_div,        Zdiv , cwhat='divergence',  iverbose=idebug )
    nPc, PDF_cnv = computePDF( xbin_bounds_div, xbin_center_div,        Zcnv , cwhat='convergence', iverbose=idebug )


    cfroot = 'PDF_'+corigin+'_dt'+cdtbin+'_'+str(reskm)+'km_'+cperiod   
    
    # Saving PDFs in `npz` files:
    np.savez_compressed( cd_in+'/'+cfroot+'_shear.npz',      name='shear',      origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nPs, xbin_bounds=xbin_bounds_shr, xbin_center=xbin_center_shr, PDF=PDF_shr )
    
    np.savez_compressed( cd_in+'/'+cfroot+'_total.npz',      name='total',      origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nPs, xbin_bounds=xbin_bounds_tot, xbin_center=xbin_center_tot, PDF=PDF_tot )
    
    np.savez_compressed( cd_in+'/'+cfroot+'_absDiv.npz', name='Divergence', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nPD, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_Div )
    
    np.savez_compressed( cd_in+'/'+cfroot+'_divergence.npz', name='divergence', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nPd, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_div )
    
    np.savez_compressed( cd_in+'/'+cfroot+'_convergence.npz', name='convergence', origin=corigin,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nPc, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_cnv )


    if l_add_gaussian:
        from math import pi, sqrt
        m_shear = np.mean(ZshrClean)
        s_shear = mjt.StdDev( m_shear, ZshrClean )
        print('LOLO: mean and sigma for shear =', m_shear, s_shear)
        zE = (xbin_center_shr[:] - m_shear )/s_shear
        #zE = xbin_center_shr[:]/s_shear

        zIntShear = np.sum( PDF_shr[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) )
        print(' Integral of PDF of shear =', zIntShear)        
        PDF_NormS = np.exp( -0.5 * zE[:]*zE[:] ) / ( s_shear * sqrt(2.*pi) )
        zIntGauss = np.sum( PDF_NormS[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) )        
        zcorr = zIntShear/zIntGauss
        print(' Integral of PDF of Gaussian =', zIntGauss)
        
        PDF_NormS = PDF_NormS*zcorr        
        print(' New Integral of PDF of Gaussian =', np.sum( PDF_NormS[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) ))

        
        np.savez_compressed( cd_in+'/'+cfroot+'_GAUSSIAN-shear.npz',      name='shear',      origin=corigin,
                             reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                             Np=nPs, xbin_bounds=xbin_bounds_shr, xbin_center=xbin_center_shr, PDF=PDF_NormS )



    
    zIntShr = np.sum( PDF_shr[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) )
    print(' * Integral of PDF of shear =', round(zIntShr,2) )        

    zIntDiv = np.sum( PDF_div[:]*(xbin_bounds_div[1:]-xbin_bounds_div[:-1]) )
    print(' * Integral of PDF of divergence =', round(zIntDiv,2),'\n')        

    zIntTot = np.sum( PDF_tot[:]*(xbin_bounds_tot[1:]-xbin_bounds_tot[:-1]) )
    print(' * Integral of PDF of tot. def. =', round(zIntTot,2) )        

        

    if iplot>0:
        cdir = './figs'
        if not path.exists(cdir): mkdir(cdir)
        
        kk = mjt.LogPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPs, name='Shear',
                            cfig=cdir+'/loglog'+cfroot+'_shear'+cxtra+'.png',      title=corigin, period=cperiod, origin=corigin )
            
        kk = mjt.LogPDFdef( xbin_bounds_tot, xbin_center_tot, PDF_tot, Np=nPs, name='Total',
                            cfig=cdir+'/loglog'+cfroot+'_total'+cxtra+'.png',      title=corigin, period=cperiod, origin=corigin )
            
        kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPD, name='|Divergence|',
                            cfig=cdir+'/loglog'+cfroot+'_Divergence'+cxtra+'.png', title=corigin, period=cperiod, origin=corigin )

        kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, Np=nPd, name='Divergence',
                            cfig=cdir+'/loglog'+cfroot+'_divergence'+cxtra+'.png', title=corigin, period=cperiod, origin=corigin )

        kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_cnv, Np=nPc, name='Convergence',
                            cfig=cdir+'/loglog'+cfroot+'_convergence'+cxtra+'.png', title=corigin, period=cperiod, origin=corigin )


        if idebug>0:
            kk = mjt.PlotPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPs, name='Shear',
                                 cfig=cdir+'/'+cfroot+'_shear'+cxtra+'.png',           title=corigin, period=cperiod )
            
            kk = mjt.PlotPDFdef( xbin_bounds_tot, xbin_center_tot, PDF_tot, Np=nPs, name='Total',
                                 cfig=cdir+'/'+cfroot+'_total'+cxtra+'.png',           title=corigin, period=cperiod )
            
            kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPD, name='Divergence',
                                 cfig=cdir+'/'+cfroot+'_Divergence'+cxtra+'.png',      title=corigin, period=cperiod )
        
            kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPd, name='Divergence',
                                 cfig=cdir+'/'+cfroot+'_divergence'+cxtra+'.png',      title=corigin, period=cperiod )
        
            kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPc, name='Convergence',
                                 cfig=cdir+'/'+cfroot+'_convergence'+cxtra+'.png',      title=corigin, period=cperiod )
    
        if l_add_gaussian:
            kk = mjt.LogPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_NormS, Np=nPs, name='Shear',
                                cfig=cdir+'/loglog'+cfroot+'_GAUSSIAN-shear'+cxtra+'.png',      title=corigin, period=cperiod )

            kk = mjt.PlotPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_NormS, Np=nPs, name='Shear',
                                 cfig=cdir+'/'+cfroot+'_GAUSSIAN_shear'+cxtra+'.png',           title=corigin, period=cperiod )
