#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, makedirs
#from glob import glob
import numpy as np
from re import split
#
import mojito   as mjt
from mojito.util import epoch2clock as e2c
from mojito import config as cfg

idebug=1
iffrmt='png'
#iffrmt='svg'
NameArcticProj='SmallArctic'

quality_mode = 'rgps'

if __name__ == '__main__':

    if not len(argv) in [3,4]:
        print('Usage: '+argv[0]+' <vardef> <def_FIELD.npz> (<def_FIELD2.npz>)')
        exit(0)

    cv_in = argv[1]
    cf_in = argv[2]
    l2files = (len(argv)==4)

    if not cv_in in ['divergence','shear','total']:
        print('ERROR: wrong deformation variable:', cv_in, '!!!'); exit(0)

    mjt.chck4f(cf_in)
    with np.load(cf_in) as data:
        FD    = data[cv_in]
        itime = int(data['time'])
        nP    = data['Npoints']
        Xc    = data['Xc']
        Yc    = data['Yc']
        X4    = data['X4']
        Y4    = data['Y4']
        corig = str(data['origin'])
        reskm = int(data['reskm_nmnl'])
    print(' * date =', e2c(itime))
    print(' * orig =', corig)
    print(' * reskm =', reskm)
    print(' * nP =', nP)

    k1 = cfg.initialize(                mode=quality_mode )
    k2 = cfg.updateConfig4Scale( reskm, mode=quality_mode )
    
    if l2files:
        cf_in2 = argv[3]        
        mjt.chck4f(cf_in2)
        with np.load(cf_in2) as data:
            FD2    = data[cv_in]
            itime2 = int(data['time'])
            nP2    = data['Npoints']
            Xc2    = data['Xc']
            Yc2    = data['Yc']
            corig2 = str(data['origin'])
            reskm2 = int(data['reskm_nmnl'])
        print(' * date2 =', e2c(itime2))
        print(' * orig2 =', corig2)
        print(' * reskm2 =', reskm2)
        print(' * nP2 =', nP2)

            
    

    # Simplifying name of experiments in certain cases:
    # NEMO-SI3_NANUK4 -> SI3
    corig = str.replace( corig, 'NEMO-SI3_NANUK4', 'SI3' )
    corig = str.replace( corig, '_', '-')
    cSclKM = str(reskm)+'km'
    if l2files:
        corig2 = str.replace( corig2, 'NEMO-SI3_NANUK4', 'SI3' )
        corig2 = str.replace( corig2, '_', '-' )
        if not reskm2==reskm: cSclKM = ''
        
    cfxtraScl, cnxtraScl, cscale = '', '', ''
    if cSclKM != '':
        cfxtraScl = '_'+cSclKM
        cscale    = 'scale = '+cSclKM
        cnxtraScl = ', '+cscale
        


        
    fdir = './figs/maps'
    if not path.exists(fdir): makedirs( fdir, exist_ok=True )


    zrx = [ np.min(Xc)-25. , np.max(Xc)+25. ]
    zry = [ np.min(Yc)-25. , np.max(Yc)+25. ]
    
    cresinfo = '('+str(reskm)+' km)'

    nmproj=NameArcticProj



    if cv_in == 'divergence':
        cwhat = 'div'
        
    elif cv_in == 'shear':
        cwhat = 'shr'
        mjt.ShowDefQuadGeoArctic( X4, Y4, cfg.rc_day2sec*FD, cfig=fdir+'/map_'+cv_in+'.'+iffrmt,
                                  nmproj=NameArcticProj, cwhat=cwhat,
                                  pFmin=0., pFmax=cfg.rc_shr_max_fig, zoom=1,
                                  rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corig+': '+cv_in+' '+cresinfo, idate=itime )
        
    elif cv_in == 'total':
        cwhat = 'tot'




    
    # Filled quads projected on the Arctic map:


    exit(0)














    
    if   l2files:

        cfroot = 'Comp_PDF_'+corig+'_vs_'+corig2+'_'+cfname+'_'+cperiod+cfxtraScl
        
        # Only log-log !
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=fdir+'/loglog'+cfroot+'.'+iffrmt, reskm=reskm,
                            title=cName+': '+cscale, period=cperiod, origin=corig,
                            ppdf2=PDF2, Np2=nP2, origin2=corig2 )    
    
    else:
        # log-log and histogram:
        cfroot = 'PDF_'+corig+'_'+cfname+'_'+cperiod+cfxtraScl
        
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=fdir+'/loglog'+cfroot+'.'+iffrmt, reskm=reskm,
                            title=cName+cnxtraScl, period=cperiod )    
    
        kk = mjt.PlotPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=fdir+'/'+cfroot+'.'+iffrmt,
                             title=cName+': '+corig+cnxtraScl, period=cperiod )
    
