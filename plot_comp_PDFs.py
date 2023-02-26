#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
from glob import glob
import numpy as np
from re import split
#
from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=1
iffrmt='png'


if __name__ == '__main__':

    if not len(argv) in [2,3]:
        print('Usage: '+argv[0]+' <PDFfile.npz> (<PDFfile2.npz>)')
        exit(0)
    cf_in = argv[1]

    l2files = (len(argv)==3)

    
    mjt.chck4f(cf_in)
    with np.load(cf_in) as data:
        cname = str(data['name'])
        corig = str(data['origin'])
        cperiod = str(data['period'])
        nP    = data['Np']
        xbin_bounds = data['xbin_bounds']
        xbin_center = data['xbin_center']
        PDF   = data['PDF']
    print(' * cname =', cname)
    print(' * cperiod =', cperiod)
    print(' * nP =', nP)


    if l2files:
        cf_in2 = argv[2]        
        mjt.chck4f(cf_in2)
        with np.load(cf_in2) as data:
            cname2 = str(data['name'])
            corig2 = str(data['origin'])
            cperiod2 = str(data['period'])
            nP2    = data['Np']
            xbin_bounds2 = data['xbin_bounds']
            xbin_center2 = data['xbin_center']
            PDF2   = data['PDF']
        print('\n * cname_2 =', cname)
        print(' * cperiod_2 =', cperiod)
        print(' * nP_2 =', nP)
        if cname2!=cname:
            print('ERROR: `cname2!=cname` !',cname2,cname)
            exit(0)
        if np.sum(np.abs(xbin_bounds2-xbin_bounds))!=0:
            print('ERROR: PDF in file 2 looks too different than in first file in terms of bin bounds?...')
            #or cperiod2!=cperiod
            exit(0)
            
    if cname == 'divergence':
        cName = 'Divergence'
    elif cname == 'shear':
        cName = 'Shear'
    elif cname == 'deftot':
        cName = 'Total deformation'
    else:
        print(' ERROR: unknow field:',cname,'!') ; exit(0)
    

    cdir = './figs'
    if not path.exists(cdir): mkdir(cdir)

    


    if l2files:

        cfroot = 'Comp_PDF_'+corig+'_vs_'+corig2+'_'+cname+'_'+cperiod
        
        # Only log-log !
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'_'+cname+'.'+iffrmt,
                            title=cName+': '+corig+' vs '+corig2, period=cperiod, origin=corig, ppdf2=PDF2, origin2=corig2 )    
    
    else:
        # log-log and histogram:
        cfroot = 'PDF_'+corig+'_'+cname+'_'+cperiod
        
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'_'+cname+'.'+iffrmt,
                            title=cName, period=cperiod )    
    
        kk = mjt.PlotPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/'+cfroot+'_'+cname+'.'+iffrmt,
                             title=cName+': '+corig, period=cperiod )
    
