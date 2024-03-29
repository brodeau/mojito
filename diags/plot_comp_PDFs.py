#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, makedirs
from glob import glob
import numpy as np
from re import split
#
import mojito   as mjt
from mojito import config as cfg

idebug=0

l_hovmuller=True

if __name__ == '__main__':

    if not len(argv) in [2,3,4,5]:
        print('Usage: '+argv[0]+' <PDFfile.npz> (<PDFfile2.npz>) (<PDFfile3.npz>) (<PDFfile3.npz>)')
        exit(0)
    cf_in = argv[1]

    l2files = (len(argv)==3)
    l3files = (len(argv)==4)
    l4files = (len(argv)==5)
    
    mjt.chck4f(cf_in)
    with np.load(cf_in) as data:
        cname = str(data['name'])
        corig = str(data['origin'])
        reskm = int(data['reskm_nmnl'])
        dtbin = int(data['dtbin'])
        #cperiod = str(data['period'])
        nP    = data['Np']
        xbin_bounds = data['xbin_bounds']
        xbin_center = data['xbin_center']
        PDF   = data['PDF']
    print(' * name =', cname)
    print(' * orig =', corig)
    print(' * reskm =', reskm)
    print(' * dtbin =', dtbin)
    #print(' * period =', cperiod)
    print(' * nP =', nP)

    if l2files or l3files or l4files:
        cf_in2 = argv[2]        
        mjt.chck4f(cf_in2)
        with np.load(cf_in2) as data:
            cname2 = str(data['name'])
            corig2 = str(data['origin'])
            reskm2 = int(data['reskm_nmnl'])
            dtbin2 = int(data['dtbin'])
            #cperiod2 = str(data['period'])
            nP2    = data['Np']
            xbin_bounds2 = data['xbin_bounds']
            xbin_center2 = data['xbin_center']
            PDF2   = data['PDF']
        print('\n * name_2 =', cname2)
        print(' * orig_2 =', corig2)
        print(' * reskm_2 =', reskm2)
        print(' * dtbin_2 =', dtbin2)
        #print(' * period_2 =', cperiod2)
        print(' * nP_2 =', nP2)
        if cname2!=cname:
            print('ERROR: `cname2!=cname` !',cname2,cname)
            exit(0)
        if np.sum(np.abs(xbin_bounds2-xbin_bounds))!=0:
            print('ERROR: PDF in file 2 looks too different than in first file in terms of bin bounds?...')
            exit(0)

    if l3files or l4files:
        cf_in3 = argv[3]        
        mjt.chck4f(cf_in3)
        with np.load(cf_in3) as data:
            cname3 = str(data['name'])
            corig3 = str(data['origin'])
            reskm3 = int(data['reskm_nmnl'])
            dtbin3 = int(data['dtbin'])
            #cperiod3 = str(data['period'])
            nP3    = data['Np']
            xbin_bounds3 = data['xbin_bounds']
            xbin_center3 = data['xbin_center']
            PDF3   = data['PDF']
        print('\n * name_3 =', cname3)
        print(' * orig_3 =', corig3)
        print(' * reskm_3 =', reskm3)
        print(' * dtbin_3 =', dtbin3)
        #print(' * period_3 =', cperiod3)
        print(' * nP_3 =', nP3)
        if cname3!=cname:
            print('ERROR: `cname3!=cname` !',cname3,cname)
            exit(0)
        if np.sum(np.abs(xbin_bounds3-xbin_bounds))!=0:
            print('ERROR: PDF in file 3 looks too different than in first file in terms of bin bounds?...')
            exit(0)

    if l4files:
        cf_in4 = argv[4]        
        mjt.chck4f(cf_in4)
        with np.load(cf_in4) as data:
            cname4 = str(data['name'])
            corig4 = str(data['origin'])
            reskm4 = int(data['reskm_nmnl'])
            dtbin4 = int(data['dtbin'])
            #cperiod4 = str(data['period'])
            nP4    = data['Np']
            xbin_bounds4 = data['xbin_bounds']
            xbin_center4 = data['xbin_center']
            PDF4   = data['PDF']
        print('\n * name_4 =', cname4)
        print(' * orig_4 =', corig4)
        print(' * reskm_4 =', reskm4)
        print(' * dtbin_4 =', dtbin4)
        #print(' * period_4 =', cperiod4)
        print(' * nP_4 =', nP4)
        if cname4!=cname:
            print('ERROR: `cname4!=cname` !',cname4,cname)
            exit(0)
        if np.sum(np.abs(xbin_bounds4-xbin_bounds))!=0:
            print('ERROR: PDF in file 4 looks too different than in first file in terms of bin bounds?...')
            exit(0)

            
            
    if   cname == 'Divergence':
        cName = '|Divergence|';      cfname = 'AbsDiv'
    elif cname == 'divergence':
        cName = 'Divergence'  ;      cfname = cName
    elif cname == 'convergence':
        cName = 'Convergence' ;      cfname = cName
    elif cname == 'shear':
        cName = 'Shear';             cfname = cName
    elif cname == 'total':
        cName = 'Total Deformation'; cfname = 'TotDef'    
    else:
        print(' ERROR: unknow field:',cname,'!') ; exit(0)
    
    

    # Simplifying name of experiments in certain cases:
    # NEMO-SI3_NANUK4 -> SI3
    corig = str.replace( corig, 'NEMO-SI3_NANUK4', 'SI3' )
    corig = str.replace( corig, '_', '-')
    cSclKM = str(reskm)+'km'
    if l2files or l3files or l4files:
        corig2 = str.replace( corig2, 'NEMO-SI3_NANUK4', 'SI3' )
        corig2 = str.replace( corig2, '_', '-' )
        if not reskm2==reskm: cSclKM = ''
    if l3files or l4files:
        corig3 = str.replace( corig3, 'NEMO-SI3_NANUK4', 'SI3' )
        corig3 = str.replace( corig3, '_', '-' )
        if not( reskm2==reskm and reskm3==reskm ): cSclKM = ''
    if l4files:
        corig4 = str.replace( corig4, 'NEMO-SI3_NANUK4', 'SI3' )
        corig4 = str.replace( corig4, '_', '-' )
        if not( reskm2==reskm and reskm4==reskm ): cSclKM = ''
        
    cfxtraScl, cnxtraScl, cscale = '', '', ''
    if cSclKM != '':
        cfxtraScl = '_'+cSclKM
        cscale    = 'scale = '+cSclKM
        cnxtraScl = ', '+cscale
        


    k2 = cfg.updateConfig4Scale( reskm, mode='rgps' )
        
    fdir = './figs/PDFs'
    if not path.exists(fdir): makedirs( fdir, exist_ok=True )


    figsfx = mjt.GetFigDotSuffix()

    
    if   l2files:

        cfroot = 'Comp_PDF_'+corig+'_vs_'+corig2+'_'+cfname+'_'+cfxtraScl
        
        # Only log-log !
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=fdir+'/loglog'+cfroot+figsfx, reskm=reskm,
                            title=cName+' ('+cscale+')', period=cperiod, origin=corig,
                            ppdf2=PDF2, Np2=nP2, origin2=corig2 )    
    
    elif l3files:

        cfroot = 'Comp_PDF_'+corig+'_vs_'+corig2+'_vs_'+corig3+'_'+cfname+'_'+cfxtraScl

        if not cfg.lc_StrictPDF:
            cfroot += '_Min1e-3'

        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, name=cName, cfig=fdir+'/loglog'+cfroot+figsfx, reskm=reskm,
                            lShowSlope=True, lAddHovm=l_hovmuller,
                            title=cName, origin=corig,
                            ppdf2=PDF2,  origin2=corig2,
                            ppdf3=PDF3,  origin3=corig3 )

        #Diag a la Bouchat:
        #kk = mjt.HovmlrPDF( xbin_bounds, xbin_center, PDF, name=cName, cfig=fdir+'/Hovmuller_PDF'+cfroot+figsfx, reskm=reskm, lShowSlope=True,
        #                    title=cName, origin=corig,
        #                    ppdf2=PDF2,  origin2=corig2,
        #                    ppdf3=PDF3,  origin3=corig3 )
        
    elif l4files:

        cfroot = 'Comp_PDF_'+corig+'_vs_'+corig2+'_vs_'+corig4+'_'+cfname+'_'+cfxtraScl

        if not cfg.lc_StrictPDF:
            cfroot += '_Min1e-3'
        
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, name=cName, cfig=fdir+'/loglog'+cfroot+figsfx, reskm=reskm,
                            title=cName, origin=corig,
                            ppdf2=PDF2,  origin2=corig2,
                            ppdf3=PDF3,  origin3=corig3,
                            ppdf4=PDF4,  origin4=corig4 )
    
    else:
        # log-log and histogram:
        cfroot = 'PDF_'+corig+'_'+cfname+'_'+cfxtraScl
        
        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=fdir+'/loglog'+cfroot+figsfx, reskm=reskm,
                            title=cName+cnxtraScl, period=cperiod )    
    
        kk = mjt.PlotPDFdef( xbin_bounds, xbin_center, PDF, Np=nP, name=cName, cfig=fdir+'/'+cfroot+figsfx,
                             title=cName+': '+corig+cnxtraScl, period=cperiod )
    
