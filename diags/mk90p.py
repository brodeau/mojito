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




def PlotP90Series( vt, VS, figname='fig_series_P90.png', y_range=(0.,0.1), x_range=None, dy=0.01 ):
    '''
        * vt: time array as Epoch Unix time 
    '''
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime as dt
    #ras = np.mean(VS) ; ram = np.mean(VM)
    #ymin, ymax, dy = sym_round_bounds(min(np.min(VM-ram),np.min(VS-ras)), max(np.max(VM-ram), np.max(VS-ras)), base=0.1 )


    if x_range:
        (xmin,xmax) = x_range
    else:
        (xmin,xmax) = (np.min(vt), np.max(vt))


    (nT,) = np.shape(vt)
    vdate = np.zeros(nT, dtype='U10')

    for jt in range(nT):
        #cdate = 
        #print(cdate)
        vdate[jt] = e2c(vt[jt], precision='D')
        #vdate = e2c(vt[jt])
        #print(vdate[jt])
    #exit(0)

    
    VX = [dt.strptime(d, "%Y-%m-%d").date() for d in vdate]


        
    fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')

    ax  = plt.axes([0.07, 0.22, 0.9, 0.75])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.xticks(rotation='60')

    plt.plot(VX, vt*0.0         , '-', color='k',                     label=None,  zorder=5)
    plt.plot(VX, VS, 'o', color='r', markersize=6, alpha=0.5, label='RGPS', zorder=10)
    #plt.plot(VX, VM-ram, '.', color=clr_mod, markersize=6, alpha=0.5, label=clabM, zorder=15)

    (ymin,ymax) = y_range


    plt.yticks( np.arange(ymin, ymax+dy, dy) )
    ax.set_ylim(ymin,ymax)
    ##ax.set_xlim(VX[0],VX[-1])
    
    #plt.ylabel('SSH [m]')
    ax.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(bbox_to_anchor=(0.55, 1.), ncol=1, shadow=True, fancybox=True)

    #plt.savefig( figname, dpi=120, transparent=False)
    print(' *** Saving figure',figname)
    plt.savefig( figname )
    plt.close(1)
    return 0





if __name__ == '__main__':

    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <file_deformation_gatheres.npz>')
        exit(0)

    cf_in = argv[1]

    # What fields are we dealing with based on filename:
    cv_in = split( '_', path.basename(cf_in) )[1]
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

    
    cd_in = path.dirname(cf_in)
    
    with np.load(cf_in) as data:
        dtbin = data['dtbin']
        reskm = data['reskm_nmnl']
        corigin = str(data['origin'])        
        cperiod = str(data['period'])
        nbF = int(data['Nbatch'])
        vdates_batch = data['dates_batch']
        Zdat = data['dates_point']
        ZDEF = data[cvar]
        

    cdtbin = str(dtbin)
    reskm  = reskm
    creskm = str(reskm)



    if np.shape(vdates_batch) != (nbF,):
        print('ERROR: problem #0!'); exit(0)


    VF90 = np.zeros( nbF )
    imsk = np.zeros( nbF, dtype='i1' )
    
    # For large scales, we must increase the size of bins:
    #zmin_div, zmin_shr, zmin_tot = 0.003, 0.003, 0.003 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
    #if   reskm >= 15. and reskm < 35.:
    #    zmin_div, zmin_shr, zmin_tot = 0.001, 0.001, 0.001
    #elif   reskm >= 35. and reskm < 50.:
    #    zmin_div, zmin_shr, zmin_tot = 5.e-4, 5.e-4, 5.e-4
    #    wVbin_min = 1.5*wVbin_min        
    #    rfexp_bin = rfexp_bin*1.2
    #elif reskm >= 50. and reskm < 100.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-4, 1.e-4, 1.e-4
    #    wVbin_min = 2.*wVbin_min
    #    rfexp_bin = rfexp_bin*1.5
    #elif reskm >= 100. and reskm < 200.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-5, 1.e-5, 1.e-5
    #    wVbin_min = 4.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.
    #elif reskm >= 200. and reskm < 400.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-5, 1.e-5, 1.e-5
    #    wVbin_min = 6.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.5
    #elif reskm >= 400. and reskm < 700.:
    #    zmin_div, zmin_shr, zmin_tot = 1.e-5, 1.e-5, 1.e-5
    #    wVbin_min = 8.*wVbin_min
    #    rfexp_bin = rfexp_bin*2.5
    #if reskm >= 35.:
    #    print('\n * Increased the width of bins and exp growth!!! Because large scale! wVbin_min, rfexp_bin =',wVbin_min, rfexp_bin)
            
    #print(ZDEF)


    print('\n *** All availabled dates for deformation:')

    ic = 0
    for jd in vdates_batch:
        print('\n *** '+e2c(jd))

        (idxDate,) = np.where( Zdat == jd )
        nV = len(idxDate)
        print('         => '+str(nV)+' '+cv_in+' deformation for this date....')
        
        ztmp = ZDEF[idxDate]
        #xdef = np.sort(ztmp[::-1])
        xdef = np.sort(ztmp)

        if xdef.shape != (nV,):
            print('ERROR: problem #1'); exit(0)
                        
        #for zd in xdef:
        #    print(zd)


        ri90 = 0.9*float(nV)
        #print(" ri90 = ", ri90, '/', nV)
        i90 = int( floor( ri90 ) )
        #print(" i90 = ", i90, '/', nV)
        rw = ri90 - float(i90)
        #print(" rw = ", rw)

        #print('lolo: nV, Nmin, i90, nV-1 ')
        if nV>=Nmin and i90<nV-2:
            # sample size must be large enough
            # and: otherwize we are too close to the end of the series, this probably a bad batch!
            v90 = (1.-rw)*xdef[i90] + rw*xdef[i90+1]
            #print(" v90 =",v90, " v90_floor=",xdef[i90], " v90_ceil=",xdef[i90+1])
            
            VF90[ic] = v90
            imsk[ic] = 1
            if v90>=Nmin:
                imsk[ic] = 1

            
        print('         => v90 =', v90, '(mask =',imsk[ic],')')
            
        ic+=1


        
    print('')


    V90 = np.ma.masked_where(imsk!=1, VF90)



    kk= PlotP90Series( vdates_batch, V90, figname='./figs/series.png', y_range=(0.,0.1), dy=0.01 )
