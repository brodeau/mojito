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

lPlotClouds = True
lPlotPDFs   = True

Nmin = 1000 ; # smallest `N` (size of the sample at a given date) required to accept a value for the 90th percentile ()


def Construct90P( ifile, vdates_batch, pdates, pdef ):
    '''
        * ifile: number of the file treated
        
    '''
    lmaskArrays = True
    (nbF,) = np.shape(vdates_batch)
    
    Z90P = np.zeros( nbF )
    zdat = np.zeros( nbF, dtype=int )
    imsk = np.zeros( nbF, dtype='i1' )

    zdates, zdef = pdates.copy(), pdef.copy()

    # Need to remove erroneous extremely small values: rc_div_min, rc_shr_min, rc_tot_min !
    from mojito import config as cfg
    kk = cfg.updateConfig4Scale( 10, mode='rgps', ltalk=False )
    print(' USING: rc_tot_min =',cfg.rc_tot_min,'to clean all data!!!')        
    (idxKeep,) = np.where(pdef>cfg.rc_tot_min)
    zdef = pdef[idxKeep]
    zdates = pdates[idxKeep]

    
    ic = 0
    for jd in vdates_batch:
        print('\n *** file #'+str(ifile)+''+e2c(jd))

        (idxDate,) = np.where( zdates == jd )
        nV = len(idxDate)
        print('         => '+str(nV)+' '+cv_in+' deformation for this date....')
                
        if nV>=Nmin:
            # sample size must be large enough
            ztmp = zdef[idxDate]
            #print('LOLO: ztmp[-5:]=',ztmp[-5:])
            #zmax = np.max(ztmp)
            #ztmp[:] = ztmp[:]/zmax ; # Normalize with highest value            
            #Z90P[ic] = np.nanpercentile(ztmp, 90, 'median_unbiased')
            Z90P[ic] = np.percentile(ztmp, 90)
            #Z90P[ic] = Z90P[ic]*zmax ; # de-normalize!
            zdat[ic] =  jd
            imsk[ic] = 1
            
        ic+=1

    if lmaskArrays:
        Z90Pm = np.ma.masked_where(imsk!=1, Z90P)
        zdatm = np.ma.masked_where(imsk!=1, zdat)
    
        zt = np.ma.MaskedArray.compressed( zdatm )
        zx = zt[1:] - zt[:-1]
        if np.any(zx<=0):
            print('ERROR [Construct90P]: `zdatm` is not increasing!'); exit(0)
        print('')        
        return zt, np.ma.MaskedArray.compressed( Z90Pm )
    else:
        (idxM,) = np.where(imsk==0)
        Z90P[idxM] = 0.

        return zdat, Z90P





def Construct90P_old( ifile, vdates_batch, pdates, pdef ):
    '''
        * ifile: number of the file treated
        
    '''
    lmaskArrays = True
    (nbF,) = np.shape(vdates_batch)
    
    Z90P = np.zeros( nbF )
    zdat = np.zeros( nbF, dtype=int )
    imsk = np.zeros( nbF, dtype='i1' )

    zdates, zdef = pdates.copy(), pdef.copy()

    ic = 0
    for jd in vdates_batch:
        print('\n *** file #'+str(ifile)+''+e2c(jd))

        (idxDate,) = np.where( zdates == jd )
        nV = len(idxDate)
        print('         => '+str(nV)+' '+cv_in+' deformation for this date....')
        
        ztmp = zdef[idxDate]
        zdfw = np.sort(ztmp) ; # Sorted in increasing order!
        #for rr in zdfw:
        #    print(rr)
        #exit(0)
        
        if zdfw.shape != (nV,):
            print('ERROR [Construct90P_old()]: problem #1'); exit(0)

        ri90 = 0.9*float(nV)
        i90f = int( floor( ri90 ) ) - 1   ; # `-1` because C-indexing...
        rw = ri90 - float(i90f)           ; # weight for interpolation

        if nV>=Nmin and i90f<nV-2:
            # sample size must be large enough
            # and: otherwize we are too close to the end of the series, this probably a bad batch!
            z90 = (1.-rw)*zdfw[i90f] + rw*zdfw[i90f+1]
            
            Z90P[ic] = z90
            zdat[ic] =  jd
            imsk[ic] = 1
            if z90>=Nmin:
                imsk[ic] = 1            
            #print('         => z90 =', z90, '(mask =',imsk[ic],')')
            
        ic+=1

    if lmaskArrays:
        Z90Pm = np.ma.masked_where(imsk!=1, Z90P)
        zdatm = np.ma.masked_where(imsk!=1, zdat)
    
        zt = np.ma.MaskedArray.compressed( zdatm )
        zx = zt[1:] - zt[:-1]
        if np.any(zx<=0):
            print('ERROR [Construct90P_old]: `zdatm` is not increasing!'); exit(0)
        print('')        
        return zt, np.ma.MaskedArray.compressed( Z90Pm )
    else:
        (idxM,) = np.where(imsk==0)
        Z90P[idxM] = 0.

        return zdat, Z90P

    

def PlotP90Series( vt1,V1, vt2=[],V2=[], vt3=[],V3=[], field=None, dt_days=3,
                   figname='fig_series_P90.png', y_range=(0.,0.1), x_range=None, dy=0.01 ):
    '''
        * vt1: time array as Epoch Unix time 
    '''
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime as dt

    vcol = mjt.vcolor
    vlw  = mjt.vlwdth
    vorg = mjt.vorig

    (nT1,) = np.shape(vt1)

    if np.shape(V1) != (nT1,):
        print('ERROR: `V1` and `vt1` do not agree in shape!'); exit(0)

    ldo2 = ( np.shape(V2)!=(0,) and np.shape(vt2)==np.shape(V2) )
    ldo3 = ( np.shape(V3)!=(0,) and np.shape(vt3)==np.shape(V3) )

    if ldo3 and not ldo2:
        print('ERROR: `ldo3 and not ldo2` !!!'); exit(0)
    
    if x_range:
        (xmin,xmax) = x_range
    else:
        (xmin,xmax) = (np.min(vt1), np.max(vt1))
        
    #vdate1 = np.array([ e2c(vt1[jt], precision='D') for jt in range(nT1) ], dtype='U10')    
    #VX1 = [dt.strptime(d, "%Y-%m-%d").date() for d in vdate1]
    vdate1 = np.array([ e2c(dt) for dt in vt1 ], dtype='U19')    
    VX1    = np.array([dt.strptime(d, "%Y-%m-%d_%H:%M:%S").date() for d in vdate1])
    
    # Construt the generic time axis for the figure => PRETTY MUCH USELESS!!!!
    dt_sec = dt_days*24*3600
    vtg = np.arange(xmin,xmax+dt_sec,dt_sec)
    vdg = np.array([ e2c(vtg[jt], precision='D') for jt in range(len(vtg)) ], dtype='U10')    
    VX  = [dt.strptime(d, "%Y-%m-%d").date() for d in vdg]
    
    if ldo2:
        (nT2,) = np.shape(vt2)
        vdate2 = np.array([ e2c(vt2[jt], precision='D') for jt in range(nT2) ], dtype='U10')    
        VX2 = [dt.strptime(d, "%Y-%m-%d").date() for d in vdate2]
        
    if ldo3:
        (nT3,) = np.shape(vt3)
        vdate3 = np.array([ e2c(vt3[jt], precision='D') for jt in range(nT3) ], dtype='U10')    
        VX3 = [dt.strptime(d, "%Y-%m-%d").date() for d in vdate3]
    
    kk = mjt.initStyle( fntzoom=1.5 )
        
    fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')
    ax  = plt.axes([0.07, 0.18, 0.9, 0.75])

    # Date ticks stuff:
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.xticks(rotation='60')
    #ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    #
    #ax.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))
    #ax.xaxis.set_minor_locator(mdates.MonthLocator())
    #ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=(1,5,10), interval=3, tz=None))
    
    plt.plot(VX1, V1, 'o-', color=vcol[0], linewidth=vlw[0], markersize=6, alpha=0.9,     label=vorg[0], zorder=10)
    if ldo2:
        plt.plot(VX2, V2, 'o-', color=vcol[1], linewidth=vlw[1], markersize=6, alpha=0.9, label=vorg[1], zorder=10)
    if ldo3:
        plt.plot(VX3, V3, 'o-', color=vcol[2], linewidth=vlw[2], markersize=6, alpha=0.9, label=vorg[2], zorder=10)

    (ymin,ymax) = y_range

    plt.yticks( np.arange(ymin, ymax+dy, dy) )
    ax.set_ylim(ymin,ymax)

    if field:
        plt.ylabel(r'P90: '+field+' [days$^{-1}$]')
    ax.grid(color='0.5', linestyle='-', linewidth=0.3)
    plt.legend(bbox_to_anchor=(0.55, 1.), ncol=1, shadow=True, fancybox=True)

    print(' *** Saving figure',figname)
    plt.savefig( figname )
    plt.close(1)
    return 0




def PlotCloud( ki, xdate, xdef, field=None, dt_days=3,
               figname='fig_series_P90.png', y_range=(0.,0.1), x_range=None, dy=0.01, zoom=1 ):
    '''
        * 
    '''
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime as dt

    vcol = mjt.vcolor
    vlw  = mjt.vlwdth
    vorg = mjt.vorig

    if np.shape(xdate) != np.shape(xdef):
        print('ERROR: `xdate` and `xdef` do not agree in shape!'); exit(0)

    if x_range:
        (xmin,xmax) = x_range
    else:
        (xmin,xmax) = (np.min(xdate), np.max(xdate))

    vdate = np.array([ e2c(dt) for dt in xdate ], dtype='U19')    
    VX    = np.array([dt.strptime(d, "%Y-%m-%d_%H:%M:%S").date() for d in vdate])
    
    kk = mjt.initStyle( fntzoom=1.5*zoom )
        
    fig = plt.figure(num = 1, figsize=(zoom*20,zoom*7), facecolor='w', edgecolor='k')
    ax  = plt.axes([0.07, 0.18, 0.9, 0.75])
    
    #plt.scatter( VX, xdef, marker=',', c=vcol[ki-1], alpha=0.3, s=1, linewidths=0, edgecolors='r' )
    plt.scatter( VX, xdef, marker='o', c=vcol[ki-1], alpha=0.2, s=2, linewidths=0, edgecolors='r' )
    #, cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, *, edgecolors=None, plotnonfinite=False, data=None, **kwargs)
    # , s=None, c=None


    # Date ticks stuff:
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.xticks(rotation='60')
    #ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())

    (ymin,ymax) = y_range
    plt.yticks( np.arange(ymin, ymax+dy, dy) )
    ax.set_ylim(ymin,ymax)

    #if field:
    #    plt.ylabel(r'P90: '+field+' [days$^{-1}$]')
    #ax.grid(color='0.5', linestyle='-', linewidth=0.3)
    #plt.legend(bbox_to_anchor=(0.55, 1.), ncol=1, shadow=True, fancybox=True)

    print(' *** Saving figure',figname)
    plt.savefig( figname, dpi=300 )
    plt.close(1)
    return 0








if __name__ == '__main__':

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
        print('ERROR: wrong `cv_in` !!! ', cv_in); exit(0)

    mjt.chck4f(cf_in1)
    with np.load(cf_in1) as data:
        dtbin1 = data['dtbin']
        reskm1 = data['reskm_nmnl']
        corigin1 = str(data['origin'])        
        #cperiod1 = str(data['period'])
        nbF1 = int(data['Nbatch'])
        VDTB1 = data['dates_batch']
        Zdat1 = data['dates_point']
        ZDEF1 = data[cvar]
    cdtbin1 = str(dtbin1)
    if np.shape(VDTB1) != (nbF1,):
        print('ERROR: problem #0 for file 1!'); exit(0)

    #if corigin1=='RGPS':
    #    # Need to remove erroneous extremely small values: rc_div_min, rc_shr_min, rc_tot_min !
    #    from mojito import config as cfg
    #    kk = cfg.updateConfig4Scale( 10, mode='rgps', ltalk=False )
    #    print(' USING: rc_tot_min =',cfg.rc_tot_min,'to clean RGPS data!!!')        
    #    (idxKeep,) = np.where(ZDEF1>cfg.rc_tot_min)
    #    ZDEF1 = ZDEF1[idxKeep]
    #    Zdat1 = Zdat1[idxKeep]

    if ldo2 or ldo3:
        cf_in2 = argv[2]
        mjt.chck4f(cf_in2)        
        with np.load(cf_in2) as data:
            dtbin2 = data['dtbin']
            reskm2 = data['reskm_nmnl']
            corigin2 = str(data['origin'])        
            #cperiod2 = str(data['period'])
            nbF2 = int(data['Nbatch'])
            VDTB2 = data['dates_batch']
            Zdat2 = data['dates_point']
            ZDEF2 = data[cvar]
            if reskm2!=reskm1:
                print('ERROR: wrong resolutin for '+cf_in2+'!!! ', reskm2); exit(0)
        cdtbin2 = str(dtbin2)
        if np.shape(VDTB2) != (nbF2,):
            print('ERROR: problem #0 for file 2!'); exit(0)
        if corigin2=='RGPS':
            print('ERROR: oops did not expect second file to be RGPS!'); exit(0)
        if np.shape(VDTB2)!=np.shape(VDTB1):
            print('ERROR: `shape(VDTB2)!=shape(VDTB1)`!'); exit(0)
        

    if ldo3:
        cf_in3 = argv[3]
        mjt.chck4f(cf_in3)        
        with np.load(cf_in3) as data:
            dtbin3 = data['dtbin']
            reskm3 = data['reskm_nmnl']
            corigin3 = str(data['origin'])        
            nbF3 = int(data['Nbatch']) # 
            VDTB3 = data['dates_batch']
            Zdat3 = data['dates_point']
            ZDEF3 = data[cvar]
            if reskm3!=reskm1:
                print('ERROR: wrong resolutin for '+cf_in3+'!!! ', reskm3); exit(0)
        cdtbin3 = str(dtbin3)
        if np.shape(VDTB3) != (nbF3,):
            print('ERROR: problem #0 for file 3!'); exit(0)
        if corigin3=='RGPS':
            print('ERROR: oops did not expect third file to be RGPS!'); exit(0)
        if np.shape(VDTB3)!=np.shape(VDTB1):
            print('ERROR: `shape(VDTB3)!=shape(VDTB1)`!'); exit(0)

        
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
                        
    if lPlotClouds:
        k0 = PlotCloud( 1, Zdat1, ZDEF1, field=cfield, figname='./figs/Cloud_'+corigin1+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )


        
    if cfield=='shear' and lPlotPDFs and ldo3:
        max_shr = 1.5 ; # day^-1
        zmin_div, zmin_shr, zmin_tot = 0.003, 0.003, 0.003 ; # day^-1
        rfexp_bin = 0.25
        wVbin_min = 0.0005 ; # Narrowest bin width (for the smalles values of deformation)

        (Nt,) = np.shape(VDTB1)
        for jt in range(Nt):
            kd1, kd2, kd3 = VDTB1[jt],VDTB2[jt],VDTB3[jt] ;  # date
            print('\n *** Preparing PDF plot batch at date:',e2c(kd1))
            (idxDate1,) = np.where( Zdat1 == kd1 )
            (idxDate2,) = np.where( Zdat2 == kd2 )
            (idxDate3,) = np.where( Zdat3 == kd3 )
            nV1,nV2,nV3 = len(idxDate1),len(idxDate2),len(idxDate3)
            print('         => '+str(nV1)+' '+cv_in+' deformation for this date....')

            if nV1>=Nmin and nV2>=Nmin and nV3>=Nmin:
                zdf1, zdf2, zdf3 = ZDEF1[idxDate1], ZDEF2[idxDate2], ZDEF3[idxDate3]

            nBinsS, xbin_bounds, xbin_center = mjt.constructExpBins( rfexp_bin, zmin_shr, max_shr, wVbin_min, name='shear' )
            nP1, PDF1 = mjt.computePDF( xbin_bounds, xbin_center,    zdf1,  cwhat='shear' )
            nP2, PDF2 = mjt.computePDF( xbin_bounds, xbin_center,    zdf2,  cwhat='shear' )
            nP3, PDF3 = mjt.computePDF( xbin_bounds, xbin_center,    zdf3,  cwhat='shear' )

            #kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF1, Np=nP1, name='shear', cfig='./figs/PDF_shear_'+e2c(kd1)+'.png' )
            #, reskm=reskm,
            #title=cName+cnxtraScl, period=cperiod )

            kk = mjt.LogPDFdef( xbin_bounds, xbin_center, PDF1, Np=nP1, name='shear', cfig='./figs/PDF_shear_'+e2c(kd1,precision='D')+'.png',
                                origin=mjt.vorig[0], title='Shear', period=e2c(kd1),
                                ppdf2=PDF2, Np2=nP2, origin2=mjt.vorig[1], ppdf3=PDF3, Np3=nP3, origin3=mjt.vorig[2] )

            



    VDAT1, V90P1 = Construct90P(1, VDTB1, Zdat1, ZDEF1 )
    
    if ldo3:
        if lPlotClouds:
            k0 = PlotCloud( 2, Zdat2, ZDEF2, field=cfield, figname='./figs/Cloud_'+corigin2+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )
            k0 = PlotCloud( 3, Zdat3, ZDEF3, field=cfield, figname='./figs/Cloud_'+corigin3+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )
        cfig = 'fig_series_P90_RGPS-BBM-aEVP_'+cfield+'.png'        
        VDAT2, V90P2 = Construct90P(2, VDTB2, Zdat2, ZDEF2 )
        VDAT3, V90P3 = Construct90P(3, VDTB3, Zdat3, ZDEF3 )

        kk= PlotP90Series( VDAT1,V90P1, vt2=VDAT2,V2=V90P2, vt3=VDAT3,V3=V90P3, field=cfield, figname='./figs/'+cfig, y_range=(0.,0.08), dy=0.01 )
        
    elif ldo2:
        if lPlotClouds:
            k0 = PlotCloud( 2, Zdat2, ZDEF2, field=cfield, figname='./figs/Cloud_'+corigin2+'_'+cfield+'.png', y_range=(0.,1.), dy=0.1, zoom=1 )
        cfig = 'fig_series_P90_RGPS-BBM'+cfield+'.png'        
        VDAT2, V90P2 = Construct90P(2, VDTB2, Zdat2, ZDEF2  )
    
    else:
        cfig = 'fig_series_P90_RGPS'+cfield+'.png'        
        kk= PlotP90Series( VDAT1, V90P1, field=cfield, figname='./figs/'+cfig, y_range=(0.,0.1), dy=0.01 )
