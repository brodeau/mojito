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


if __name__ == '__main__':

    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <PDFfile.npz>')
        exit(0)
    cf_in = argv[1]


    mjt.chck4f(cf_in)

    with np.load(cf_in) as data:
        cname = str(data['name'])
        cperiod = str(data['period'])
        wbin    = data['wbin']
        nP    = data['Np']
        xbin_bounds = data['xbin_bounds']
        xbin_center = data['xbin_center']
        PDF   = data['PDF']

    print(' * cname =', cname)
    print(' * cperiod =', cperiod)
    print(' * nP =', nP)



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




    

        
    kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/loglog'+cfroot+'_'+cname+'.svg',
                        wbin=wbin, title=cname, period=cperiod )
    
    xfrng=[0.001,0.1] ; # x-range we want on the x-axis of the plot
    
    kk = mjt.PlotPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=cdir+'/'+cfroot+'_'+cname+'.svg',
                         xrng=xfrng, wbin=wbin, title=cname, period=cperiod )
    
