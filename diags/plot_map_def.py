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
izoom = 1.5

NameArcticProj='SmallArctic'

quality_mode = 'rgps_map'

if __name__ == '__main__':

    if not len(argv) in [3,4,5]:
        print('Usage: '+argv[0]+' <vardef> <def_FIELD.npz> (<def_FIELD2.npz>) (<def_FIELD3.npz>)')
        exit(0)

    cv_in = argv[1]
    cf_in = argv[2]
    l2files = (len(argv)==4 or len(argv)==5)
    l3files = (len(argv)==5)
    
    if not cv_in in ['divergence','shear','total']:
        print('ERROR: wrong deformation variable:', cv_in, '!!!'); exit(0)

    ik = cfg.controlModeName( path.basename(__file__), quality_mode )
    k1 = cfg.initialize( mode=quality_mode )
        
    mjt.chck4f(cf_in)
    with np.load(cf_in) as data:
        FD    = cfg.rc_day2sec*data[cv_in]
        itime = int(data['time'])
        nP    = data['Npoints']
        Xc    = data['Xc']
        Yc    = data['Yc']
        X4p    = data['X4']
        Y4p    = data['Y4']
        corig = str(data['origin'])
        reskm = int(data['reskm_nmnl'])
    print('\n * date =', e2c(itime))
    print(' * orig =', corig)
    print(' * reskm =', reskm)
    print(' * nP =', nP)

    k2 = cfg.updateConfig4Scale( reskm, mode=quality_mode )
    
    if l2files:
        cf_in2 = argv[3]        
        mjt.chck4f(cf_in2)
        with np.load(cf_in2) as data:
            FD2    = cfg.rc_day2sec*data[cv_in]
            itime2 = int(data['time'])
            nP2    = data['Npoints']
            #Xc2    = data['Xc']
            #Yc2    = data['Yc']
            X4p2   = data['X4']
            Y4p2   = data['Y4']            
            corig2 = str(data['origin'])
            reskm2 = int(data['reskm_nmnl'])
        print('\n * date2 =', e2c(itime2))
        print(' * orig2 =', corig2)
        print(' * reskm2 =', reskm2)
        print(' * nP2 =', nP2)

    if l3files:
        cf_in3 = argv[4]        
        mjt.chck4f(cf_in3)
        with np.load(cf_in3) as data:
            FD3    = cfg.rc_day2sec*data[cv_in]
            itime3 = int(data['time'])
            nP3    = data['Npoints']
            #Xc3    = data['Xc']
            #Yc3    = data['Yc']
            X4p3   = data['X4']
            Y4p3   = data['Y4']                        
            corig3 = str(data['origin'])
            reskm3 = int(data['reskm_nmnl'])
        print('\n * date3 =', e2c(itime3))
        print(' * orig3 =', corig3)
        print(' * reskm3 =', reskm3)
        print(' * nP3 =', nP3)

            

    # Simplifying name of experiments in certain cases:
    # NEMO-SI3_NANUK4 -> SI3
    corig = str.replace( corig, 'NEMO-SI3_NANUK4', 'SI3' )
    corig = str.replace( corig, '_', '-')
    cSclKM = str(reskm)+'km'
    if l2files:
        corig2 = str.replace( corig2, 'NEMO-SI3_NANUK4', 'SI3' )
        corig2 = str.replace( corig2, '_', '-' )
        if not reskm2==reskm: cSclKM = ''
    if l3files:
        corig3 = str.replace( corig3, 'NEMO-SI3_NANUK4', 'SI3' )
        corig3 = str.replace( corig3, '_', '-' )
        if not reskm3==reskm: cSclKM = ''
        
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
        if not l2files:
            mjt.ShowDefQuadGeoArctic( X4p, Y4p, FD, cfig=fdir+'/map_'+cv_in+'_'+corig+'.'+iffrmt,
                                      nmproj=NameArcticProj, cwhat=cwhat,
                                      pFmin=-cfg.rc_div_max_fig, pFmax=cfg.rc_div_max_fig, zoom=izoom,
                                      rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                      title=mjt.vorig[0]+': '+cv_in+' '+cresinfo, idate=itime )
        elif l3files:
            mjt.ShowMultiDefQuadGeoArctic( X4p, Y4p, FD, X4p2, Y4p2, FD2, X4p3, Y4p3, FD3, zoom=izoom,
                                           cfig=fdir+'/map_'+cv_in+'_'+corig+'.'+iffrmt,
                                           nmproj=NameArcticProj, cwhat=cwhat,
                                           pFmin=-cfg.rc_div_max_fig, pFmax=cfg.rc_tot_max_fig,
                                           rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                           title1=mjt.vorig[0]+': '+cv_in+' '+cresinfo, idate=itime,
                                           title2=mjt.vorig[1]+': '+cv_in+' '+cresinfo, title3=mjt.vorig[2]+': '+cv_in+' '+cresinfo )
        else:
            print('Fixme!')
        
    elif cv_in == 'shear':
        cwhat = 'shr'
        if not l2files:
            mjt.ShowDefQuadGeoArctic( X4p, Y4p, FD, cfig=fdir+'/map_'+cv_in+'_'+corig+'.'+iffrmt,
                                      nmproj=NameArcticProj, cwhat=cwhat,
                                      pFmin=0., pFmax=cfg.rc_shr_max_fig, zoom=izoom,
                                      rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                      title=mjt.vorig[0]+': '+cv_in+' '+cresinfo, idate=itime )

        elif l3files:
            mjt.ShowMultiDefQuadGeoArctic( X4p, Y4p, FD, X4p2, Y4p2, FD2, X4p3, Y4p3, FD3, zoom=izoom,
                                           cfig=fdir+'/map_'+cv_in+'_'+corig+'.'+iffrmt,
                                           nmproj=NameArcticProj, cwhat=cwhat,
                                           pFmin=0., pFmax=cfg.rc_tot_max_fig,
                                           rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                           title1=mjt.vorig[0]+': '+cv_in+' '+cresinfo, idate=itime,
                                           title2=mjt.vorig[1]+': '+cv_in+' '+cresinfo, title3=mjt.vorig[2]+': '+cv_in+' '+cresinfo )
        else:
            print('Fixme!')
        
    elif cv_in == 'total':
        cwhat = 'tot'
        if not l2files:
            mjt.ShowDefQuadGeoArctic( X4p, Y4p, FD, cfig=fdir+'/map_'+cv_in+'_'+corig+'.'+iffrmt,
                                      nmproj=NameArcticProj, cwhat=cwhat,
                                      pFmin=0., pFmax=cfg.rc_tot_max_fig, zoom=izoom,
                                      rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                      title=corig+': '+cv_in+' '+cresinfo, idate=itime )
        elif l3files:
            mjt.ShowMultiDefQuadGeoArctic( X4p, Y4p, FD, X4p2, Y4p2, FD2, X4p3, Y4p3, FD3, zoom=izoom,
                                           cfig=fdir+'/map_'+cv_in+'_'+corig+'.'+iffrmt,
                                           nmproj=NameArcticProj, cwhat=cwhat,
                                           pFmin=0., pFmax=cfg.rc_tot_max_fig,
                                           rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                           title1=corig+': '+cv_in+' '+cresinfo, idate=itime,
                                           title2=mjt.vorig[1]+': '+cv_in+' '+cresinfo, title3=mjt.vorig[2]+': '+cv_in+' '+cresinfo )

        
        else:
            print('Fixme!')

 
