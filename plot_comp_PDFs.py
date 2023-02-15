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


if __name__ == '__main__':

    if not len(argv) in [2,3]:
        print('Usage: '+argv[0]+' <PDFfile.npz> (<PDFfile2.npz>)')
        exit(0)
    cf_in = argv[1]

    l2files = (len(argv)==3)

    
    mjt.chck4f(cf_in)
    with np.load(cf_in) as data:
        cname = str(data['name'])
        corigin = str(data['origin'])
        cperiod = str(data['period'])
        wbin    = data['wbin']
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
            corigin2 = str(data['origin'])
            cperiod2 = str(data['period'])
            wbin2    = data['wbin']
            nP2    = data['Np']
            xbin_bounds2 = data['xbin_bounds']
            xbin_center2 = data['xbin_center']
            PDF2   = data['PDF']
        print('\n * cname_2 =', cname)
        print(' * cperiod_2 =', cperiod)
        print(' * nP_2 =', nP)
        if cname2!=cname or np.sum(np.abs(xbin_bounds2-xbin_bounds))!=0:
            print('ERROR: PDF in file 2 looks too different than in first file...')
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

    cfroot = 'PDF_'+cname+'_'+cperiod



    xfrng=[0.001,1.] ; # x-range we want on the x-axis of the plot



    if l2files:
        # Only log-log !
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'_'+cname+'.svg',
                            wbin=wbin, title=cName+': '+corigin+' vs '+corigin2, period=cperiod, origin=corigin, ppdf2=PDF2, origin2=corigin2 )    
    

    else:
        # log-log and histogram:
        
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'_'+cname+'.svg',
                            wbin=wbin, title=cName, period=cperiod )    
    
        kk = mjt.PlotPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/'+cfroot+'_'+cname+'.svg',
                             xrng=xfrng, wbin=wbin, title=cName+': '+corigin, period=cperiod )
    
