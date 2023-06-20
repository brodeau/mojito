#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split
import mojito   as mjt
from mojito.util import epoch2clock as e2c
#
from math import floor

idebug=1
iplot=1
izoom=3

lPlotClouds = False
lPlotPDFs   = False

Nmin = 1000 ; # smallest `N` (size of the sample at a given date) required to accept a value for the 90th percentile ()

zP=98

iffrmt='png'

if __name__ == '__main__':

    ldo1 = True
    
    Narg = len(argv)
    if not Narg in [2,3,4]:
        print('Usage: '+argv[0]+' <file_deformation_gatherd.npz> (<file_deformation_gatherd.npz>) (<file_deformation_gatherd.npz>)')
        exit(0)

    cf_in1 = argv[1]

    ldo2 = (Narg==3)
    ldo3 = (Narg==4)

    # What fields are we dealing with based on filename:
    cv_in = split( '_', path.basename(cf_in1) )[1]
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
        mjt.printEE('wrong `cv_in` !!! ', cv_in)

    mjt.chck4f(cf_in1)
    with np.load(cf_in1) as data:
        dtbin1 = data['dtbin']
        dates_batch1 = data['dates_batch']
        reskm1 = data['reskm_nmnl']
        corigin1 = str(data['origin'])        
        #cperiod1 = str(data['period'])
        nbF1 = int(data['Nbatch'])
        VDTB1 = data['dates_batch']
        Zdat1 = data['dates_point']
        ZDEF1 = data[cvar]
    cdtbin1 = str(dtbin1)
    if np.shape(VDTB1) != (nbF1,):
        mjt.printEE('problem #0 for file 1!')

    if ldo2 or ldo3:
        cf_in2 = argv[2]
        mjt.chck4f(cf_in2)        
        with np.load(cf_in2) as data:
            dtbin2 = data['dtbin']
            dates_batch2 = data['dates_batch']
            reskm2 = data['reskm_nmnl']
            corigin2 = str(data['origin'])        
            #cperiod2 = str(data['period'])
            nbF2 = int(data['Nbatch'])
            VDTB2 = data['dates_batch']
            Zdat2 = data['dates_point']
            ZDEF2 = data[cvar]
            if reskm2!=reskm1:
                mjt.printEE('wrong resolutin for '+cf_in2+'!!! ', reskm2)
        cdtbin2 = str(dtbin2)
        if np.shape(VDTB2) != (nbF2,):
            mjt.printEE('problem #0 for file 2!')
        if corigin2=='RGPS':
            mjt.printEE('oops did not expect second file to be RGPS!')
        if np.shape(VDTB2)!=np.shape(VDTB1):
            print('    ', np.shape(VDTB2), np.shape(VDTB1))
            print(dates_batch2)
            print(dates_batch1)
            mjt.printW('`shape(VDTB2)!=shape(VDTB1)`!')
        

    if ldo3:
        cf_in3 = argv[3]
        mjt.chck4f(cf_in3)        
        with np.load(cf_in3) as data:
            dtbin3 = data['dtbin']
            dates_batch3 = data['dates_batch']
            reskm3 = data['reskm_nmnl']
            corigin3 = str(data['origin'])        
            nbF3 = int(data['Nbatch']) # 
            VDTB3 = data['dates_batch']
            Zdat3 = data['dates_point']
            ZDEF3 = data[cvar]
            if reskm3!=reskm1:
                mjt.printEE('wrong resolutin for '+cf_in3+'!!! ', reskm3)
        cdtbin3 = str(dtbin3)
        if np.shape(VDTB3) != (nbF3,):
            mjt.printEE('problem #0 for file 3!')
        if corigin3=='RGPS':
            mjt.printEE('oops did not expect third file to be RGPS!')
        if np.shape(VDTB3)!=np.shape(VDTB1):
            print('    ', np.shape(VDTB3), np.shape(VDTB1))
            mjt.printW('`shape(VDTB3)!=shape(VDTB1)`!')

        
    reskm  = reskm1
    creskm = str(reskm)


    if cfield=='divergence':
        idxD  = np.where(ZDEF1>0.)
        ZDEF1 = ZDEF1[idxD]
        Zdat1 = Zdat1[idxD]
        #
        if ldo2 or ldo3:
            idxD  = np.where(ZDEF2>0.)
            ZDEF2 = ZDEF2[idxD]
            Zdat2 = Zdat2[idxD]
        if ldo3:
            idxD  = np.where(ZDEF3>0.)
            ZDEF3 = ZDEF3[idxD]
            Zdat3 = Zdat3[idxD]


    if ldo3:
        print(' * Number of points for '+mjt.vorig[0]+', '+mjt.vorig[1]+', '+mjt.vorig[2]+' =',len(ZDEF1),len(ZDEF2),len(ZDEF3))


            
    if lPlotClouds:
        k0 = mjt.PlotCloud( 1, Zdat1, ZDEF1, field=cfield, figname='./figs/Cloud_'+corigin1+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )


    ##DEBUG:
    #if cfield=='shear' and lPlotPDFs and ldo3:
    #    max_shr = 1.5 ; # day^-1
    #    zmin_div, zmin_shr, zmin_tot = 0.003, 0.003, 0.003 ; # day^-1
    #    rfexp_bin = 0.25
    #    wVbin_min = 0.0005 ; # Narrowest bin width (for the smalles values of deformation)
    #
    #    (Nt,) = np.shape(VDTB1)
    #    for jt in range(Nt):
    #        kd1, kd2, kd3 = VDTB1[jt],VDTB2[jt],VDTB3[jt] ;  # date
    #        print('\n *** Preparing PDF plot batch at date:',e2c(kd1))
    #        (idxDate1,) = np.where( Zdat1 == kd1 )
    #        (idxDate2,) = np.where( Zdat2 == kd2 )
    #        (idxDate3,) = np.where( Zdat3 == kd3 )
    #        nV1,nV2,nV3 = len(idxDate1),len(idxDate2),len(idxDate3)
    #        print('         => '+str(nV1)+' '+cv_in+' deformation for this date....')
    #
    #        if nV1>=Nmin and nV2>=Nmin and nV3>=Nmin:
    #            zdf1, zdf2, zdf3 = ZDEF1[idxDate1], ZDEF2[idxDate2], ZDEF3[idxDate3]
    #
    #        nBinsS, xbin_bounds, xbin_center = mjt.constructExpBins( rfexp_bin, zmin_shr, max_shr, wVbin_min, name='shear' )
    #        nP1, PDF1 = mjt.computePDF( xbin_bounds, xbin_center,    zdf1,  cwhat='shear' )
    #        nP2, PDF2 = mjt.computePDF( xbin_bounds, xbin_center,    zdf2,  cwhat='shear' )
    #        nP3, PDF3 = mjt.computePDF( xbin_bounds, xbin_center,    zdf3,  cwhat='shear' )
    #
    #        #kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF1, Np=nP1, name='shear', cfig='./figs/PDF_shear_'+e2c(kd1)+'.png' )
    #        #, reskm=reskm,
    #        #title=cName+cnxtraScl, period=cperiod )
    #
    #        kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF1, Np=nP1, name='shear', cfig='./figs/PDF_shear_'+e2c(kd1,precision='D')+'.png',
    #                            origin=mjt.vorig[0], title='Shear', period=e2c(kd1),
    #                            ppdf2=PDF2, Np2=nP2, origin2=mjt.vorig[1], ppdf3=PDF3, Np3=nP3, origin3=mjt.vorig[2] )
    ##DEBUG.
            



    VDAT1, V90P1 = mjt.Construct90P(1, VDTB1, Zdat1, ZDEF1, pp=zP, Nmin=Nmin )
    
    if ldo3:

        if lPlotClouds:
            k0 = mjt.PlotCloud( 2, Zdat2, ZDEF2, field=cfield, figname='./figs/Cloud_'+corigin2+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )
            k0 = mjt.PlotCloud( 3, Zdat3, ZDEF3, field=cfield, figname='./figs/Cloud_'+corigin3+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )
        cfig = 'fig_series_P'+str(zP)+'_RGPS-BBM-aEVP_'+cfield+'.'+iffrmt        
        VDAT2, V90P2 = mjt.Construct90P(2, VDTB2, Zdat2, ZDEF2, pp=zP, Nmin=Nmin )
        VDAT3, V90P3 = mjt.Construct90P(3, VDTB3, Zdat3, ZDEF3, pp=zP, Nmin=Nmin )

        kk= mjt.PlotP90Series( VDAT1,V90P1, vt2=VDAT2,V2=V90P2, vt3=VDAT3,V3=V90P3, field=cfield, whatP=str(zP), figname='./figs/'+cfig, y_range=(0.,0.25), dy=0.05 )
        
    elif ldo2:
        if lPlotClouds:
            k0 = mjt.PlotCloud( 2, Zdat2, ZDEF2, field=cfield, figname='./figs/Cloud_'+corigin2+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )
        cfig = 'fig_series_P'+str(zP)+'_RGPS-BBM'+cfield+'.'+iffrmt        
        VDAT2, V90P2 = mjt.Construct90P(2, VDTB2, Zdat2, ZDEF2, pp=zP, Nmin=Nmin  )
    
    else:
        cfig = 'fig_series_P'+str(zP)+'_RGPS'+cfield+'.'+iffrmt        
        kk= mjt.PlotP90Series( VDAT1, V90P1, field=cfield, figname='./figs/'+cfig, y_range=(0.,0.1), dy=0.01 )
