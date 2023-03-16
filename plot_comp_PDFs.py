#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
from glob import glob
import numpy as np
from re import split
#
import mojito   as mjt

idebug=1
#iffrmt='png'
iffrmt='svg'


if __name__ == '__main__':

    if not len(argv) in [2,3,4]:
        print('Usage: '+argv[0]+' <PDFfile.npz> (<PDFfile2.npz>) (<PDFfile3.npz>)')
        exit(0)
    cf_in = argv[1]

    l2files = (len(argv)==3)
    l3files = (len(argv)==4)
    
    mjt.chck4f(cf_in)
    with np.load(cf_in) as data:
        cname = str(data['name'])
        corig = str(data['origin'])
        reskm = int(data['reskm_nmnl'])
        cperiod = str(data['period'])
        nP    = data['Np']
        xbin_bounds = data['xbin_bounds']
        xbin_center = data['xbin_center']
        PDF   = data['PDF']
    print(' * name =', cname)
    print(' * orig =', corig)
    print(' * reskm =', reskm)
    print(' * period =', cperiod)
    print(' * nP =', nP)

    if l2files or l3files:
        cf_in2 = argv[2]        
        mjt.chck4f(cf_in2)
        with np.load(cf_in2) as data:
            cname2 = str(data['name'])
            corig2 = str(data['origin'])
            reskm2 = int(data['reskm_nmnl'])
            cperiod2 = str(data['period'])
            nP2    = data['Np']
            xbin_bounds2 = data['xbin_bounds']
            xbin_center2 = data['xbin_center']
            PDF2   = data['PDF']
        print('\n * name_2 =', cname)
        print(' * orig_2 =', corig)
        print(' * reskm_2 =', reskm)        
        print(' * period_2 =', cperiod)
        print(' * nP_2 =', nP)
        if cname2!=cname:
            print('ERROR: `cname2!=cname` !',cname2,cname)
            exit(0)
        if np.sum(np.abs(xbin_bounds2-xbin_bounds))!=0:
            print('ERROR: PDF in file 2 looks too different than in first file in terms of bin bounds?...')
            #or cperiod2!=cperiod
            exit(0)

    if l3files:
        cf_in3 = argv[3]        
        mjt.chck4f(cf_in3)
        with np.load(cf_in3) as data:
            cname3 = str(data['name'])
            corig3 = str(data['origin'])
            reskm3 = int(data['reskm_nmnl'])
            cperiod3 = str(data['period'])
            nP3    = data['Np']
            xbin_bounds3 = data['xbin_bounds']
            xbin_center3 = data['xbin_center']
            PDF3   = data['PDF']
        print('\n * name_3 =', cname)
        print(' * orig_3 =', corig)
        print(' * reskm_3 =', reskm)
        print(' * period_3 =', cperiod)
        print(' * nP_3 =', nP)
        if cname3!=cname:
            print('ERROR: `cname3!=cname` !',cname3,cname)
            exit(0)
        if np.sum(np.abs(xbin_bounds3-xbin_bounds))!=0:
            print('ERROR: PDF in file 3 looks too different than in first file in terms of bin bounds?...')
            #or cperiod3!=cperiod
            exit(0)

    exit(0)
            
        
    if   cname == 'Divergence':
        cName = '|Divergence|'
    elif cname == 'divergence':
        cName = 'Divergence'
    elif cname == 'convergence':
        cName = 'Convergence'
    elif cname == 'shear':
        cName = 'Shear'
    elif cname == 'deftot':
        cName = 'Total deformation'
    else:
        print(' ERROR: unknow field:',cname,'!') ; exit(0)
    
    cfname = cName
    if cname == 'Divergence': cfname = 'AbsDiv'    
        
    cdir = './figs'
    if not path.exists(cdir): mkdir(cdir)

    
    if   l2files:

        cfroot = 'Comp_PDF_'+corig+'_vs_'+corig2+'_'+cfname+'_'+cperiod
        
        # Only log-log !
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'.'+iffrmt,
                            title=cName+': '+corig+' vs '+corig2, period=cperiod, origin=corig,
                            ppdf2=PDF2, Np2=nP2, origin2=corig2 )    
    
    elif l3files:

        cfroot = 'Comp_PDF_'+corig+'_vs_'+corig2+'_vs_'+corig3+'_'+cfname+'_'+cperiod
        
        # Only log-log !
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'.'+iffrmt,
                            title=cName+': '+corig+' vs '+corig2+' vs '+corig2, period=cperiod, origin=corig,
                            ppdf2=PDF2, Np2=nP2, origin2=corig2, ppdf3=PDF3, Np3=nP3, origin3=corig3 )
    
    else:
        # log-log and histogram:
        cfroot = 'PDF_'+corig+'_'+cfname+'_'+cperiod
        
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'.'+iffrmt,
                            title=cName, period=cperiod )    
    
        kk = mjt.PlotPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/'+cfroot+'.'+iffrmt,
                             title=cName+': '+corig, period=cperiod )
    
