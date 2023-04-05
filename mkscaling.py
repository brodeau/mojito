#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir, environ
#from glob import glob
import numpy as np
#from re import split
import mojito   as mjt

idebug=1
iplot=1

dir_in = '<HOME>/Nextcloud/data/mojitoNEW'

l_cst_bins = False ; rfexp_bin = 0.2

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

do_scales = [ 10, 20, 40, 320, 640 ]



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










    
    exit(0)

    
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
        
    # Control agreement:
    if dtbin1!=dtbin2:
        print('ERROR: `dtbin1!=dtbin2` !!!'); exit(0)
    if reskm1!=reskm2:
        print('ERROR: `reskm1!=reskm2` !!!'); exit(0)
    if corig1!=corig2:
        print('ERROR: `corig1!=corig2` !!!'); exit(0)
    if cperd1!=cperd2:
        print('ERROR: `cperd1!=cperd2` !!!'); exit(0)
    if nbF1!=nbF2:
        print('ERROR: `nbF1!=nbF2` !!!'); exit(0)

    dtbin  = dtbin1
    cdtbin = str(dtbin1)
    reskm  = reskm1
    creskm = str(reskm1)
    corigin = corig1
    cperiod = cperd1

    # For large scales, we must increase the size of bins:
    #min_div = 0.003 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
    #min_shr = 0.003 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
    #if   reskm >= 15. and reskm < 35.:
    #    min_div, min_shr = 0.001, 0.001
    #elif   reskm >= 35. and reskm < 50.:
    #    min_div, min_shr = 5.e-4, 5.e-4
    #    wVbin_min = 1.5*wVbin_min        
    #    rfexp_bin = rfexp_bin*1.2
    #elif reskm >= 50. and reskm < 100.:
    #    min_div, min_shr = 1.e-4, 1.e-4
    #    wVbin_min = 2.*wVbin_min
    #    rfexp_bin = rfexp_bin*1.5
    #elif reskm >= 100. and reskm < 200.:
    #    min_div, min_shr = 1.e-5, 1.e-5
    #    wVbin_min = 4.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.
    #elif reskm >= 200. and reskm < 400.:
    #    min_div, min_shr = 1.e-5, 1.e-5
    #    wVbin_min = 6.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.5
    #elif reskm >= 400. and reskm < 700.:
    #    min_div, min_shr = 1.e-5, 1.e-5
    #    wVbin_min = 8.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.5
    #if reskm >= 35.:
    #    print('\n * Increased the width of bins and exp growth!!! Because large scale! wVbin_min, rfexp_bin =',wVbin_min, rfexp_bin)
    #        
    #print(ZDiv)

    cfroot = 'PDF_'+corigin+'_dt'+cdtbin+'_'+str(reskm)+'km_'+cperiod   
    
