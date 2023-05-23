#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
#from glob import glob
import numpy as np
from re import split
import mojito   as mjt
from mojito.util import epoch2clock as e2c
#
from math import floor

idebug=1
iplot=1
izoom=3

Nmin = 1000 ; # smallest `N` (size of the sample at a given date) required to accept a value for the 90th percentile ()




def Construct90P( ifile, vdates_batch, pdates, pdef ):
    '''
        * ifile: number of the file treated
        
    '''
    lmaskArrays = True
    (nbF,) = np.shape(vdates_batch)
    
    z90p = np.zeros( nbF )
    zdat = np.zeros( nbF, dtype=int )
    imsk = np.zeros( nbF, dtype='i1' )
    
    ic = 0
    for jd in vdates_batch:
        print('\n *** file #'+str(ifile)+''+e2c(jd))

        (idxDate,) = np.where( pdates == jd )
        nV = len(idxDate)
        print('         => '+str(nV)+' '+cv_in+' deformation for this date....')
        
        ztmp = pdef[idxDate]
        xdef = np.sort(ztmp)

        if xdef.shape != (nV,):
            print('ERROR [Construct90P()]: problem #1'); exit(0)

        ri90 = 0.9*float(nV)
        i90 = int( floor( ri90 ) )
        rw = ri90 - float(i90)

        if nV>=Nmin and i90<nV-2:
            # sample size must be large enough
            # and: otherwize we are too close to the end of the series, this probably a bad batch!
            z90 = (1.-rw)*xdef[i90] + rw*xdef[i90+1]
            
            z90p[ic] = z90
            zdat[ic] =  jd
            imsk[ic] = 1
            if z90>=Nmin:
                imsk[ic] = 1            
            #print('         => z90 =', z90, '(mask =',imsk[ic],')')
            
        ic+=1

    if lmaskArrays:
        z90pm = np.ma.masked_where(imsk!=1, z90p)
        zdatm = np.ma.masked_where(imsk!=1, zdat)
    
        zt = np.ma.MaskedArray.compressed( zdatm )
        zx = zt[1:] - zt[:-1]
        if np.any(zx<=0):
            print('ERROR [Construct90P]: `zdatm` is not increasing!'); exit(0)
        #print(zx)    
        print('')        
        return zt, np.ma.MaskedArray.compressed( z90pm )
    else:
        (idxM,) = np.where(imsk==0)
        z90p[idxM] = 0.

        return zdat, z90p

    

def PlotP90Series( vt1,V1, vt2=[],V2=[], vt3=[],V3=[], field=None, figname='fig_series_P90.png', vlabels=[], y_range=(0.,0.1), x_range=None, dy=0.01 ):
    '''
        * vt1: time array as Epoch Unix time 
    '''
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime as dt
    #ras = np.mean(V1) ; ram = np.mean(VM)
    #ymin, ymax, dy = sym_round_bounds(min(np.min(VM-ram),np.min(V1-ras)), max(np.max(VM-ram), np.max(V1-ras)), base=0.1 )

    import mojito.plot as mjtp

    vcol = mjtp.vcolor
    vlw  = mjtp.vlwdth
    vorg = mjtp.vorig

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

    vdate1 = np.array([ e2c(vt1[jt], precision='D') for jt in range(nT1) ], dtype='U10')    
    VX1 = [dt.strptime(d, "%Y-%m-%d").date() for d in vdate1]

    if ldo2:
        (nT2,) = np.shape(vt2)
        vdate2 = np.array([ e2c(vt2[jt], precision='D') for jt in range(nT2) ], dtype='U10')    
        VX2 = [dt.strptime(d, "%Y-%m-%d").date() for d in vdate2]
        
    if ldo3:
        (nT3,) = np.shape(vt3)
        vdate3 = np.array([ e2c(vt3[jt], precision='D') for jt in range(nT3) ], dtype='U10')    
        VX3 = [dt.strptime(d, "%Y-%m-%d").date() for d in vdate3]
        
    
    fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')
    ax  = plt.axes([0.07, 0.18, 0.9, 0.75])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.xticks(rotation='60')

    plt.plot(VX1, vt1*0.0         , '-', color='k',                     label=None,  zorder=5)
    plt.plot(VX1, V1, 'o-', color=vcol[0], linewidth=vlw[0], markersize=6, alpha=0.9,     label=vorg[0], zorder=10)
    if ldo2:
        plt.plot(VX2, V2, 'o-', color=vcol[1], linewidth=vlw[1], markersize=6, alpha=0.9, label=vorg[1], zorder=10)
    if ldo3:
        plt.plot(VX3, V3, 'o-', color=vcol[2], linewidth=vlw[2], markersize=6, alpha=0.9, label=vorg[2], zorder=10)

    (ymin,ymax) = y_range


    plt.yticks( np.arange(ymin, ymax+dy, dy) )
    ax.set_ylim(ymin,ymax)
    ##ax.set_xlim(VX1[0],VX1[-1])

    if field:
        plt.ylabel(r'P90: '+field+' [days$^{-1}$]')
    ax.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(bbox_to_anchor=(0.55, 1.), ncol=1, shadow=True, fancybox=True)

    #plt.savefig( figname, dpi=120, transparent=False)
    print(' *** Saving figure',figname)
    plt.savefig( figname )
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

    
    #cd_in = path.dirname(cf_in1)
    
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
                
    if ldo2 or ldo3:
        cf_in2 = argv[2]
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

    if ldo3:
        cf_in3 = argv[3]
        with np.load(cf_in3) as data:
            dtbin3 = data['dtbin']
            reskm3 = data['reskm_nmnl']
            corigin3 = str(data['origin'])        
            #cperiod3 = str(data['period'])
            nbF3 = int(data['Nbatch'])
            VDTB3 = data['dates_batch']
            Zdat3 = data['dates_point']
            ZDEF3 = data[cvar]
            if reskm3!=reskm1:
                print('ERROR: wrong resolutin for '+cf_in3+'!!! ', reskm3); exit(0)
        cdtbin3 = str(dtbin3)
        if np.shape(VDTB3) != (nbF3,):
            print('ERROR: problem #0 for file 3!'); exit(0)

        
    reskm  = reskm1
    creskm = str(reskm)
    

    VDAT1, V90P1 = Construct90P(1, VDTB1, Zdat1, ZDEF1  )



    


    
    if ldo3:
        cfig = 'fig_series_P90_RGPS-BBM-aEVP.png'        
        VDAT2, V90P2 = Construct90P(2, VDTB2, Zdat2, ZDEF2 )
        VDAT3, V90P3 = Construct90P(3, VDTB3, Zdat3, ZDEF3 )

        kk= PlotP90Series( VDAT1,V90P1, vt2=VDAT2,V2=V90P2, vt3=VDAT3,V3=V90P3, field=cfield, figname='./figs/'+cfig, y_range=(0.,0.08), dy=0.01 )
        
    elif ldo2:
        cfig = 'fig_series_P90_RGPS-BBM.png'        
        VDAT2, V90P2 = Construct90P(2, VDTB2, Zdat2, ZDEF2  )
    
    else:
        cfig = 'fig_series_P90_RGPS.png'        
        kk= PlotP90Series( VDAT1, V90P1, figname='./figs/'+cfig, y_range=(0.,0.1), dy=0.01 )
