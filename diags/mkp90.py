#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
#from glob import glob
import numpy as np
from re import split
import mojito   as mjt
from mojito.util import epoch2clock as e2c
#
from math import floor

idebug=1
iplot=1
izoom=3

Nmin = 500 ; # smallest `N` (size of the sample at a given date) required to accept a value for the 90th percentile ()

#l_cst_bins = False ; rfexp_bin = 0.3
l_cst_bins = False ; rfexp_bin = 0.25

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

    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <file_deformation_gatheres.npz>')
        exit(0)

    cf_in = argv[1]

    # What fields are we dealing with based on filename:
    cv_in = split( '_', path.basename(cf_in) )[1]
    print(cv_in)

    if cv_in == 'DIV':
        cfield = 'divergence'
        cvar   = 'xdiv'
    elif cv_in == 'SHR':
        cfield = 'shear'
        cvar   = 'xshr'
    elif cv_in == 'TOT':
        cfield = 'total'
        cvar   = 'xtot'
    else:
        print('ERROR: wrong `cv_in` !!! ', cv_in); exit(0)

    
    cd_in = path.dirname(cf_in)
    
    with np.load(cf_in) as data:
        dtbin = data['dtbin']
        reskm = data['reskm_nmnl']
        corigin = str(data['origin'])        
        cperiod = str(data['period'])
        nbF = int(data['Nbatch'])
        vdates_batch = data['dates_batch']
        Zdat = data['dates_point']
        ZDEF = data[cvar]
        

    cdtbin = str(dtbin)
    reskm  = reskm
    creskm = str(reskm)



    if np.shape(vdates_batch) != (nbF,):
        print('ERROR: problem #0!'); exit(0)


    VF90 = np.zeros( nbF )
    imsk = np.zeros( nbF, dtype='i1' )
    
    # For large scales, we must increase the size of bins:
    #zmin_div, zmin_shr, zmin_tot = 0.003, 0.003, 0.003 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
    #if   reskm >= 15. and reskm < 35.:
    #    zmin_div, zmin_shr, zmin_tot = 0.001, 0.001, 0.001
    #elif   reskm >= 35. and reskm < 50.:
    #    zmin_div, zmin_shr, zmin_tot = 5.e-4, 5.e-4, 5.e-4
    #    wVbin_min = 1.5*wVbin_min        
    #    rfexp_bin = rfexp_bin*1.2
    #elif reskm >= 50. and reskm < 100.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-4, 1.e-4, 1.e-4
    #    wVbin_min = 2.*wVbin_min
    #    rfexp_bin = rfexp_bin*1.5
    #elif reskm >= 100. and reskm < 200.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-5, 1.e-5, 1.e-5
    #    wVbin_min = 4.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.
    #elif reskm >= 200. and reskm < 400.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-5, 1.e-5, 1.e-5
    #    wVbin_min = 6.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.5
    #elif reskm >= 400. and reskm < 700.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-5, 1.e-5, 1.e-5
    #    wVbin_min = 8.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.5
    #if reskm >= 35.:
    #    print('\n * Increased the width of bins and exp growth!!! Because large scale! wVbin_min, rfexp_bin =',wVbin_min, rfexp_bin)
            
    #print(ZDEF)


    print('\n *** All availabled dates for deformation:')

    ic = 0
    for jd in vdates_batch:
        print('\n *** '+e2c(jd))

        (idxDate,) = np.where( Zdat == jd )
        nV = len(idxDate)
        print('         => '+str(nV)+' '+cv_in+' deformation for this date....')
        
        ztmp = ZDEF[idxDate]
        #xdef = np.sort(ztmp[::-1])
        xdef = np.sort(ztmp)

        if xdef.shape != (nV,):
            print('ERROR: problem #1'); exit(0)
                        
        #for zd in xdef:
        #    print(zd)


        ri90 = 0.9*float(nV)
        #print(" ri90 = ", ri90, '/', nV)
        i90 = int( floor( ri90 ) )
        #print(" i90 = ", i90, '/', nV)
        rw = ri90 - float(i90)
        #print(" rw = ", rw)

        #print('lolo: nV, Nmin, i90, nV-1 ')
        if nV>=Nmin and i90<nV-2:
            # sample size must be large enough
            # and: otherwize we are too close to the end of the series, this probably a bad batch!
            v90 = (1.-rw)*xdef[i90] + rw*xdef[i90+1]
            #print(" v90 =",v90, " v90_floor=",xdef[i90], " v90_ceil=",xdef[i90+1])
            
            VF90[ic] = v90
            imsk[ic] = 1
            if v90>=Nmin:
                imsk[ic] = 1

            
        print('         => v90 =', v90, '(mask =',imsk[ic],')')
            
        ic+=1


        
    print('')





    exit(0)
    cxtra = ''
    if l_cst_bins:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = constructCstBins( zmin_div, max_div, wbin_div, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = constructCstBins( zmin_shr, max_shr, wbin_shr, name='shear',      iverbose=idebug )        
        # For the total deformation:
        nBinsS, xbin_bounds_tot, xbin_center_tot = constructCstBins( zmin_tot, max_tot, wbin_tot, name='shear',      iverbose=idebug )
        
    else:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = constructExpBins( rfexp_bin, zmin_div, max_div, wVbin_min, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = constructExpBins( rfexp_bin, zmin_shr, max_shr, wVbin_min, name='shear',      iverbose=idebug )
        # For the total deformation:
        nBinsS, xbin_bounds_tot, xbin_center_tot = constructExpBins( rfexp_bin, zmin_tot, max_tot, wVbin_min, name='total',      iverbose=idebug )
        cxtra = '_incB'



    # Signed divergence:
    (idxN,) = np.where(ZDEF<-zmin_div)
    nPn  = len(idxN)
    Zcnv = np.zeros(nPn)
    Zcnv[:] = - ZDEF[idxN] ; # We want it to be positive for the graph...
    
    (idxP,) = np.where(ZDEF> zmin_div)
    nPp  = len(idxP)
    Zdiv = np.zeros(nPp)
    Zdiv[:] = ZDEF[idxP]


    if l_add_gaussian:
        nPs, PDF_shr, ZshrClean = computePDF( xbin_bounds_shr, xbin_center_shr, Zshr, cwhat='shear', return_cleaned=True   ,  iverbose=idebug )
    else:
        nPs, PDF_shr = computePDF( xbin_bounds_shr, xbin_center_shr,    Zshr,  cwhat='shear',      iverbose=idebug )

    nPs, PDF_tot = computePDF( xbin_bounds_tot, xbin_center_tot,        Ztot,  cwhat='total',       iverbose=idebug )    
    nPD, PDF_Div = computePDF( xbin_bounds_div, xbin_center_div, np.abs(ZDEF), cwhat='Divergence',  iverbose=idebug )
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
