#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
#from glob import glob
import numpy as np
 #from re import split
import mojito   as mjt
from mojito import config as cfg

idebug=1
iplot=1

#l_cst_bins = False ; rfexp_bin = 0.3
l_cst_bins = False ; rfexp_bin = 0.25

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

# About width of bins:
if l_cst_bins:
    wbin_div = 0.001 ; # day^-1
    wbin_shr = 0.001 ; # day^-1
    wbin_tot = 0.001 ; # day^-1
else:
    wVbin_min = 0.0005 ; # Narrowest bin width (for the smalles values of deformation)

l_add_gaussian = False


if __name__ == '__main__':

    #if not len(argv) in [3,5]:
    #    print('Usage: '+argv[0]+' <directory_input_npz_files> <dtbin_h> <creskm> <string_id_origin>')
    #    print('   or: '+argv[0]+' <directory_input_npz_files> <file_prefix>')
    #    exit(0)
    if not len(argv) in [3]:
        print('Usage: '+argv[0]+' <file_divergence.npz> <mode (rgps,model,xlose)>')
        exit(0)

    cf_div_in    = argv[1]

    quality_mode = argv[2]
    k1 = cfg.controlModeName( path.basename(__file__), quality_mode )

    cf_shr_in = str.replace( cf_div_in, 'DIV', 'SHR' )
    cf_tot_in = str.replace( cf_div_in, 'DIV', 'TOT' )

    cd_in = path.dirname(cf_div_in)

    with np.load(cf_div_in) as data:
        dtbin1 = data['dtbin']
        reskm1 = data['reskm_nmnl']
        corig1 = str(data['origin'])
        #cperd1 = str(data['period'])
        nbF1 = int(data['Nbatch'])
        ZDiv = data['xdiv']

    with np.load(cf_shr_in) as data:
        dtbin2 = data['dtbin']
        reskm2 = data['reskm_nmnl']
        corig2 = str(data['origin'])
        #cperd2 = str(data['period'])
        nbF2   = int(data['Nbatch'])
        Zshr   =     data['xshr']

    with np.load(cf_tot_in) as data:
        dtbin3 = data['dtbin']
        reskm3 = data['reskm_nmnl']
        corig3 = str(data['origin'])
        #cperd3 = str(data['period'])
        nbF3   = int(data['Nbatch'])
        Ztot   =     data['xtot']


    # Control agreement:
    if dtbin1!=dtbin2 or dtbin1!=dtbin3:
        print('ERROR: `dtbin1!=dtbin2` !!!'); exit(0)
    if reskm1!=reskm2 or reskm1!=reskm3:
        print('ERROR: `reskm1!=reskm2` !!!'); exit(0)
    if corig1!=corig2 or corig1!=corig3:
        print('ERROR: `corig1!=corig2` !!!'); exit(0)
    #if cperd1!=cperd2 or cperd1!=cperd3:
    #    print('ERROR: `cperd1!=cperd2` !!!'); exit(0)
    if nbF1!=nbF2 or nbF1!=nbF3:
        print('ERROR: `nbF1!=nbF2` !!!'); exit(0)

    dtbin  = dtbin1
    cdtbin = str(dtbin1)
    reskm  = reskm1
    creskm = str(reskm1)
    corigin = corig1
    #cperiod = cperd1

    # Minimum and maximum deformation values acceptable for pdf:
    k2 = cfg.updateConfig4Scale( reskm, mode='rgps', ltalk=False )
    print(' * [Construct90P()] => uses `rc_MinDef` =',cfg.rc_MinDef,'to clean all data / rgps, regardless of origin!!!')
    zmin_def, zmax_def = cfg.rc_def_min_pdf, cfg.rc_def_max_pdf  ; # day^-1

    # For large scales, we must increase the size of bins:
    if   reskm >= 35. and reskm < 50.:
        wVbin_min = 1.5*wVbin_min
        rfexp_bin = rfexp_bin*1.2
    elif reskm >= 50. and reskm < 100.:
        wVbin_min = 2.*wVbin_min
        rfexp_bin = rfexp_bin*1.5
    elif reskm >= 100. and reskm < 200.:
        wVbin_min = 4.*wVbin_min
        rfexp_bin = rfexp_bin*2.
    elif reskm >= 200. and reskm < 400.:
        wVbin_min = 6.*wVbin_min
        rfexp_bin = rfexp_bin*2.5
    elif reskm >= 400. and reskm < 700.:
        wVbin_min = 8.*wVbin_min
        rfexp_bin = rfexp_bin*2.5
    if reskm >= 35.:
        print('\n * Increased the width of bins and exp growth!!! Because large scale! wVbin_min, rfexp_bin =',wVbin_min, rfexp_bin)

    print(ZDiv)


    cxtra = ''
    if l_cst_bins:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = mjt.constructCstBins( zmin_def, zmax_def, wbin_div, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = mjt.constructCstBins( zmin_def, zmax_def, wbin_shr, name='shear',      iverbose=idebug )
        # For the total deformation:
        nBinsS, xbin_bounds_tot, xbin_center_tot = mjt.constructCstBins( zmin_def, zmax_def, wbin_tot, name='shear',      iverbose=idebug )

    else:
        # For the divergence
        nBinsD, xbin_bounds_div, xbin_center_div = mjt.constructExpBins( rfexp_bin, zmin_def, zmax_def, wVbin_min, name='divergence', iverbose=idebug )
        # For the shear:
        nBinsS, xbin_bounds_shr, xbin_center_shr = mjt.constructExpBins( rfexp_bin, zmin_def, zmax_def, wVbin_min, name='shear',      iverbose=idebug )
        # For the total deformation:
        nBinsS, xbin_bounds_tot, xbin_center_tot = mjt.constructExpBins( rfexp_bin, zmin_def, zmax_def, wVbin_min, name='total',      iverbose=idebug )
        cxtra = '_incB'



    # Signed divergence:
    (idxN,) = np.where(ZDiv<-zmin_def)
    nPn  = len(idxN)
    Zcnv = np.zeros(nPn)
    Zcnv[:] = - ZDiv[idxN] ; # We want it to be positive for the graph...

    (idxP,) = np.where(ZDiv> zmin_def)
    nPp  = len(idxP)
    Zdiv = np.zeros(nPp)
    Zdiv[:] = ZDiv[idxP]


    if l_add_gaussian:
        nPs, PDF_shr, ZshrClean = mjt.computePDF( xbin_bounds_shr, xbin_center_shr, Zshr, cwhat='shear', return_cleaned=True   ,  iverbose=idebug )
    else:
        nPs, PDF_shr = mjt.computePDF( xbin_bounds_shr, xbin_center_shr,    Zshr,  cwhat='shear',      iverbose=idebug )

    nPs, PDF_tot = mjt.computePDF( xbin_bounds_tot, xbin_center_tot,        Ztot,  cwhat='total',       iverbose=idebug )
    nPD, PDF_Div = mjt.computePDF( xbin_bounds_div, xbin_center_div, np.abs(ZDiv), cwhat='Divergence',  iverbose=idebug )
    nPd, PDF_div = mjt.computePDF( xbin_bounds_div, xbin_center_div,        Zdiv , cwhat='divergence',  iverbose=idebug )
    nPc, PDF_cnv = mjt.computePDF( xbin_bounds_div, xbin_center_div,        Zcnv , cwhat='convergence', iverbose=idebug )


    #cfroot = 'PDF_'+corigin+'_dt'+cdtbin+'_'+str(reskm)+'km_'+cperiod
    cfroot = 'PDF_'+corigin+'_dt'+cdtbin+'_'+str(reskm)+'km'

    # Saving PDFs in `npz` files:
    np.savez_compressed( cd_in+'/'+cfroot+'_shear.npz',      name='shear',      origin=corigin,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nPs, xbin_bounds=xbin_bounds_shr, xbin_center=xbin_center_shr, PDF=PDF_shr )

    np.savez_compressed( cd_in+'/'+cfroot+'_total.npz',      name='total',      origin=corigin,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nPs, xbin_bounds=xbin_bounds_tot, xbin_center=xbin_center_tot, PDF=PDF_tot )

    np.savez_compressed( cd_in+'/'+cfroot+'_absDiv.npz', name='Divergence', origin=corigin,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nPD, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_Div )

    np.savez_compressed( cd_in+'/'+cfroot+'_divergence.npz', name='divergence', origin=corigin,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nPd, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_div )

    np.savez_compressed( cd_in+'/'+cfroot+'_convergence.npz', name='convergence', origin=corigin,
                         reskm_nmnl=reskm, dtbin=dtbin,
                         Np=nPc, xbin_bounds=xbin_bounds_div, xbin_center=xbin_center_div, PDF=PDF_cnv )


    if l_add_gaussian:
        from math import pi, sqrt
        m_shear = np.mean(ZshrClean)
        s_shear = mjt.StdDev( m_shear, ZshrClean )
        print(' *** mean and sigma for shear =', m_shear, s_shear)
        zE = (xbin_center_shr[:] - m_shear )/s_shear

        zIntShear = np.sum( PDF_shr[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) )
        print(' Integral of PDF of shear =', zIntShear)
        PDF_NormS = np.exp( -0.5 * zE[:]*zE[:] ) / ( s_shear * sqrt(2.*pi) )
        zIntGauss = np.sum( PDF_NormS[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) )
        zcorr = zIntShear/zIntGauss
        print(' Integral of PDF of Gaussian =', zIntGauss)

        PDF_NormS = PDF_NormS*zcorr
        print(' New Integral of PDF of Gaussian =', np.sum( PDF_NormS[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) ))


        np.savez_compressed( cd_in+'/'+cfroot+'_GAUSSIAN-shear.npz',      name='shear',      origin=corigin,
                             reskm_nmnl=reskm, dtbin=dtbin,
                             Np=nPs, xbin_bounds=xbin_bounds_shr, xbin_center=xbin_center_shr, PDF=PDF_NormS )



    zIntShr = np.sum( PDF_shr[:]*(xbin_bounds_shr[1:]-xbin_bounds_shr[:-1]) )
    print('\n * Integral of PDF of shear =', round(zIntShr,2) )
    zIntDiv = np.sum( PDF_div[:]*(xbin_bounds_div[1:]-xbin_bounds_div[:-1]) )
    print(' * Integral of PDF of divergence =', round(zIntDiv,2) )
    zIntTot = np.sum( PDF_tot[:]*(xbin_bounds_tot[1:]-xbin_bounds_tot[:-1]) )
    print(' * Integral of PDF of tot. def. =', round(zIntTot,2),'\n' )


    if iplot>0:
        cdir = './figs'
        if not path.exists(cdir): mkdir(cdir)

        kk = mjt.LogPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPs, name='Shear',
                            cfig=cdir+'/loglog'+cfroot+'_shear'+cxtra+'.png',      title=corigin, origin=corigin )

        kk = mjt.LogPDFdef( xbin_bounds_tot, xbin_center_tot, PDF_tot, Np=nPs, name='Total',
                            cfig=cdir+'/loglog'+cfroot+'_total'+cxtra+'.png',      title=corigin, origin=corigin )

        kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPD, name='|Divergence|',
                            cfig=cdir+'/loglog'+cfroot+'_Divergence'+cxtra+'.png', title=corigin, origin=corigin )

        kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, Np=nPd, name='Divergence',
                            cfig=cdir+'/loglog'+cfroot+'_divergence'+cxtra+'.png', title=corigin, origin=corigin )

        kk = mjt.LogPDFdef( xbin_bounds_div, xbin_center_div, PDF_cnv, Np=nPc, name='Convergence',
                            cfig=cdir+'/loglog'+cfroot+'_convergence'+cxtra+'.png', title=corigin, origin=corigin )


        if idebug>0:
            kk = mjt.PlotPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nPs, name='Shear',
                                 cfig=cdir+'/'+cfroot+'_shear'+cxtra+'.png',           title=corigin )

            kk = mjt.PlotPDFdef( xbin_bounds_tot, xbin_center_tot, PDF_tot, Np=nPs, name='Total',
                                 cfig=cdir+'/'+cfroot+'_total'+cxtra+'.png',           title=corigin )

            kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPD, name='Divergence',
                                 cfig=cdir+'/'+cfroot+'_Divergence'+cxtra+'.png',      title=corigin )

            kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPd, name='Divergence',
                                 cfig=cdir+'/'+cfroot+'_divergence'+cxtra+'.png',      title=corigin )

            kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_Div, Np=nPc, name='Convergence',
                                 cfig=cdir+'/'+cfroot+'_convergence'+cxtra+'.png',      title=corigin )

        if l_add_gaussian:
            kk = mjt.LogPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_NormS, Np=nPs, name='Shear',
                                cfig=cdir+'/loglog'+cfroot+'_GAUSSIAN-shear'+cxtra+'.png',      title=corigin )

            kk = mjt.PlotPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_NormS, Np=nPs, name='Shear',
                                 cfig=cdir+'/'+cfroot+'_GAUSSIAN_shear'+cxtra+'.png',           title=corigin )
