#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J xHORRORx
#SBATCH -o out_RGPS_selection_%J.out
#SBATCH -e err_RGPS_selection_%J.err
#SBATCH --time=02:55:00
#SBATCH --account=python
#SBATCH --mem=8000
##################################################################

# TO DO:
#   * Even with a Time bin of 3 days, we want the first date (for the result) to be at <YEAR>/01/01 00:00:00
#     and not <YEAR>/01/02 12:00:00 !!!
#     => intiate the search right before the early bin boundary!
#
#

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np

from re import split
from scipy import interpolate

from netCDF4 import Dataset

from climporn import epoch2clock, clock2epoch
import mojito as mjt


Nrec_min_days = 3 ; # retain only buoys with a record length >= Nrec_min_days

iplot = 1 ; # show the result in figures?
l_drop_doublons = False
rd_tol = 5. # tolerance distance in km to conclude it's the same buoy

cdt_pattern = 'YYYY-MM-DD_00:00:00' ; # pattern for dates

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!


#idebug = 2 ; vIDdbg = [ 18 , 62 , 98, 1200, 2800, 3000, 8800, 10290, 16805 ] ; l_show_IDs_fig = True ; # Debug: pick a tiny selection of IDs to follow...
#idebug = 1 ; vIDdbg = range(0,36000,100)
idebug = 2

list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

interp_1d = 0 ; # Time interpolation to fixed time axis: 0 => linear / 1 => akima

MinDistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]







if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    for cd in ['npz','figs']:
        if not path.exists('./'+cd): mkdir('./'+cd)
    if not path.exists('./figs/SELECTION'): mkdir('./figs/SELECTION')

    ####################################################################################################
    narg = len(argv)
    if not narg in [5]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YYYYMMDD1> <YYYYMMDD2> <dt_binning (hours)>')
        exit(0)
    cf_in    =     argv[1]
    cYmmd1   =     argv[2]
    cYmmd2   =     argv[3]
    idtbin_h = int(argv[4])
    ####################################################################################################
    
    dt_bin =   float(idtbin_h*3600) ; # bin width for time scanning in [s], aka time increment while scanning for valid time intervals

    cY1,  cY2  = cYmmd1[0:4], cYmmd2[0:4]
    cmm1, cmm2 = cYmmd1[4:6], cYmmd2[4:6]
    cdd1, cdd2 = cYmmd1[6:8], cYmmd2[6:8]

    cdt1 = str.replace(cdt_pattern,'YYYY',cY1) ; cdt1 = str.replace(cdt1,'MM',cmm1) ; cdt1 = str.replace(cdt1,'DD',cdd1)
    cdt2 = str.replace(cdt_pattern,'YYYY',cY2) ; cdt2 = str.replace(cdt2,'MM',cmm2) ; cdt2 = str.replace(cdt2,'DD',cdd2)

    cf_out = 'RGPS_ice_drift_'+split('_', cdt1)[0]+'_'+split('_', cdt2)[0]+'_lb.nc' ;# netCDF file to generate


    print('\n *** Date range to restrain data to:')
    print(' ==> '+cdt1+' to '+cdt2 )
    #print('      with buoys that do not pass '+cdtI+' canceled!')


    
    rdt1, rdt2 = clock2epoch(cdt1), clock2epoch(cdt2)
    print( '   ===> in epoch time: ', rdt1, 'to', rdt2 )
    print( '       ====> double check: ', epoch2clock(rdt1), 'to',  epoch2clock(rdt2))

    rdtI = rdt1 + 3600*24*Nrec_min_days ; # minimum date for a buoy to reach
    cdtI = epoch2clock(rdtI)
    print('\n *** We shall not select buoys that do not make it to at least',cdtI)
    
    # Load `distance to coast` data:
    vlon_dist, vlat_dist, xdist = mjt.LoadDist2CoastNC( fdist2coast_nc )


    # Important we want the bins to be centered on the specified dates, so:
    rdt1 = rdt1 - 0.5*dt_bin
    rdt2 = rdt2 #- 0.5*dt_bin

    
    # Build scan time axis willingly at relative high frequency (dt_bin << dt_buoy_Nmnl)
    Nt, vTbin, cTbin =   mjt.TimeBins4Scanning( rdt1, rdt2, dt_bin, iverbose=0 )
    if idebug>1:
        for jt in range(Nt):
            print('   --- jt: '+str(jt)+' => ',epoch2clock(vTbin[jt,0]),epoch2clock(vTbin[jt,1]),epoch2clock(vTbin[jt,2]))
    
    # Open, inspect the input file and load raw data:
    Np0, vtime0, vx0, vy0, vlon0, vlat0, vBIDs0 = mjt.InspectLoadData( cf_in, list_expected_var )


    if idebug>0:
        rlat_min, rlat_max = np.min(vlat0), np.max(vlat0)
        print('\n *** Southernmost & Northernmost latitudes in the dataset =>', rlat_min, rlat_max)
    
    #if idebug>0:
    #    vIDs = np.array(vIDdbg, dtype=int)
    #    cbla = '[DEBUG] IDs =>'
    #    for ii in vIDs: cbla = cbla+' '+str(ii)
    #else:
    
    vIDs, idxUNQid = np.unique(vBIDs0, return_index=True )

    Nb   = len(vIDs)
    print('\n *** We IDed '+str(Nb)+' (unique) buoys in this file...')

    
    # First we are going to get rid of all buoys that are not ... our defined period
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # To mask IDs we cancel:
    vmask     = np.zeros(Nb, dtype='i1') + 1

    # For each buoy, getting date of earliest and latest snapshot
    print('\n *** Scanning to study precise date range of this file...')

    for jb in range(Nb):
        if jb%1000==0: print('    .... '+str(jb)+' / '+str(Nb)+'...')

        jid = vIDs[jb]

        idx = np.where( vBIDs0==jid )
        vt  = vtime0[idx]
        t1, t2  = np.min(vt), np.max(vt)

        #idx_keep = np.where((vt>=rdt1) & (vt<=rdtI))
        
        
        # Get rid of buoys that pop up for the 1st time after rdt1:
        if t1>rdt1:
            if idebug>1: print('   --- excluding buoy #'+str(jb)+' with ID: ',jid,' (pops up after '+epoch2clock(rdt1)+'!)')
            vmask[jb] = 0
        if t2<rdtI:
            if idebug>1: print('   --- excluding buoy #'+str(jb)+' with ID: ',jid,' (vanishes before '+cdtI+'!)')
            vmask[jb] = 0

        if idebug>0 and Nb<=20: print('      * buoy #'+str(jb)+' (id='+str(jid)+': '+epoch2clock(t1)+' ==> '+epoch2clock(t2))


    NbR = np.sum(vmask)
    print('\n *** Only '+str(NbR)+' buoys out of '+str(Nb)+' are remaining (time range)')

    # Masking all destroyed buoys:
    vIDs      =  np.ma.masked_where( vmask==0, vIDs      )


    # Dropping canceled (masked) buoys and updating number of valid buoys:
    vIDs = np.ma.MaskedArray.compressed(vIDs)
    Nb   = len(vIDs)
    print('\n *** UPDATE: based on date selection, there are now '+str(Nb)+' buoys to follow!')



    # Final arrays have 2 dimmensions Nb & Nt (buoy ID & time record)
    xmsk = np.zeros((Nt,Nb), dtype='i1') + 1
    xlon = np.zeros((Nt,Nb))
    xlat = np.zeros((Nt,Nb))
    xX   = np.zeros((Nt,Nb))
    xY   = np.zeros((Nt,Nb))
    
    ic = -1
    for jid in vIDs:
        ic = ic + 1

        idx  = np.where( vBIDs0==jid )
        vt   = vtime0[idx] ;            # that's the time axis of this particular buoy [epoch time]
        Np   = len(vt)     ;            # mind that Np varies from 1 buoy to another...

        if np.any( vt[1:Np]-vt[0:Np-1] < 0.):
            print('ERROR: time not increasing for buoy ID:'+str(jid)+' !!!')
            exit(0)

        Nt_b = Nt

        # Does this buoy survive beyond the fixed time range or not??? lilo
        rtend = np.max(vt)
        l_shorten = ( rtend < rdt2 )
        if l_shorten:
            if idebug>1: print('    * buoy with ID: '+str(jid)+' does not make it to '+epoch2clock(rdt2), '\n      => dies at '+epoch2clock(rtend))
            for Nt_b in range(Nt):
                if vTbin[Nt_b,0] > rtend: break
            if idebug>1: print('      ==> last index of vTbin to use is: '+str(Nt_b-1)+' => '+epoch2clock(vTbin[Nt_b-1,0]))

        if   interp_1d==0:
            fG = interpolate.interp1d(           vt, vlat0[idx])
            fC = interpolate.interp1d(           vt,   vy0[idx])
        elif interp_1d==1:
            fG = interpolate.Akima1DInterpolator(vt, vlat0[idx])
            fC = interpolate.Akima1DInterpolator(vt,   vy0[idx])
        xlat[0:Nt_b,ic] = fG(vTbin[0:Nt_b,0])
        xY[0:Nt_b,ic]   = fC(vTbin[0:Nt_b,0])
        
        if   interp_1d==0:
            fG = interpolate.interp1d(           vt, vlon0[idx])
            fC = interpolate.interp1d(           vt,   vx0[idx])
        elif interp_1d==1:
            fG = interpolate.Akima1DInterpolator(vt, vlon0[idx])
            fC = interpolate.Akima1DInterpolator(vt,   vx0[idx])
        xlon[0:Nt_b,ic] = fG(vTbin[0:Nt_b,0])
        xX[0:Nt_b,ic]   = fC(vTbin[0:Nt_b,0])
        
        if l_shorten: xmsk[Nt_b:Nt,ic] = 0

        if idebug>1 and Nb<=20:
            # Visual control of interpolation:
            lbr.plot_interp_series( jid, 'lat', vt, vTbin[0:Nt_b,0], vlat0[idx], xlat[0:Nt_b,ic] )
            lbr.plot_interp_series( jid, 'lon', vt, vTbin[0:Nt_b,0], vlon0[idx], xlon[0:Nt_b,ic] )


    #
    xlat = np.ma.masked_where( xmsk==0, xlat )
    xlon = np.ma.masked_where( xmsk==0, xlon )
    xX   = np.ma.masked_where( xmsk==0, xX   )
    xY   = np.ma.masked_where( xmsk==0, xY )

    
    # Removing doublons:
    if l_drop_doublons:
        # Use x,y !!! not lon,lat likebelow !!!
        from gonzag import Haversine

        vid_mask = np.zeros(Nb, dtype='i1')
        vid_mask[:] = 1
        
        # Based on first time step:
        jt   = 0
        vlat = xlat[jt,:]
        vlon = xlon[jt,:]
        #
        for jb in range(Nb):        
            # All the following work should be done if present buoy has not been cancelled yet
            if vid_mask[jb] == 1:
                
                jid = vIDs[jb]
                if idebug>0: print('\n *** "TOO CLOSE" scanning: Buoy #'+str(jb)+'=> ID ='+str(jid))
        
                rlat = vlat[jb]
                rlon = vlon[jb]
        
                if idebug>0: print('    ==> lat,lon = ',rlat,rlon)
                
                # Building array of distances (km) with all other buoy at this same time record:
                vdist = Haversine( rlat, rlon, vlat, vlon )
        
                # Need to mask itself (distance = 0!)
                if vdist[jb] != 0.:
                    print(' PROBLEM: distance with yourself should be 0!')
                vdist[jb] = 9999.
        
                if idebug>1:
                    for rd in vdist:
                        if rd<10.: print('   distance = ', rd)
        
                rdmin = np.min(vdist)
                if idebug>0: print('    ==> closest other buoy is at '+str(round(rdmin,2))+' km !')
                
                if rdmin < rd_tol:
                    # There is at least 1 buoy too close!
                    #  => populating the buoys closer than `rd_tol` km:
                    (vidx_close,) = np.where( vdist < rd_tol )
        
                    if idebug>0:
                        print('    ==> list of buoys closer than '+str(rd_tol)+' km:')
                        for jb in vidx_close:
                            print(vIDs[ii])
    
                    vid_mask[(vidx_close,)] = 0

        
        # Alright, all IDs spotted via vid_mask should be cancelled in the 2D fields, at all time steps
        (idmsk,) = np.where(vid_mask==0)
        xmsk[:,idmsk] = 0
        vIDs = np.ma.masked_where( vid_mask==0, vIDs )
        xlat = np.ma.masked_where( xmsk==0, xlat )
        xlon = np.ma.masked_where( xmsk==0, xlon )
        
        Nok = np.sum(xmsk[0,:])
        print('\n *** UPDATE: after "too close" canceling, there are '+str(Nok)+' buoys left at t0!')

        # Need to compress arrays, didn't find an elegant way to do it so here is the following abomination:
        xmsk2 = np.zeros((Nt,Nok), dtype=int)
        xlon2 = np.zeros((Nt,Nok))
        xlat2 = np.zeros((Nt,Nok))
        ic = 0
        for jb in range(Nb):
            if xmsk[0,jb] == 1:
                xmsk2[:,ic] = xmsk[:,jb]
                xlon2[:,ic] = xlon[:,jb]
                xlat2[:,ic] = xlat[:,jb]
                ic = ic+1
        vIDs = np.ma.MaskedArray.compressed(vIDs)
        Nb = Nok
        #print('Shape of vIDs, xlat2, xlon2, xmsk2 =',vIDs.shape,xlat2.shape,xlon2.shape,xmsk2.shape)
        del xmsk, xlat, xlon
        xmsk = xmsk2
        xlat = xlat2
        xlon = xlon2
        del xmsk2, xlat2, xlon2
        #print('Shape of vIDs, xlat, xlon, xmsk =',vIDs.shape,xlat.shape,xlon.shape,xmsk.shape)
        #exit(0)
        
    ### if l_drop_doublons

    
    xY   = np.ma.masked_where( xmsk==0, xY )
    xX   = np.ma.masked_where( xmsk==0, xX )
    xlat = np.ma.masked_where( xmsk==0, xlat )
    xlon = np.ma.masked_where( xmsk==0, xlon )
    

    # GENERATION OF COMPREHENSIVE NETCDF FILE:
    kk = mjt.ncSaveCloudBoys( cf_out, vTbin[:,0], vIDs, xY, xX, xlat, xlon, tunits=ctunits_expected )

    
    if iplot>0:

        for jt in range(Nt):
            rt = vTbin[jt,0]
            ct = epoch2clock(rt)
            print(' Ploting for '+ct+'!')
            cfig = './figs/SELECTION/t-interp_buoys_RGPS_'+ct+'.png'
            mjt.ShowBuoysMap( vTbin[jt,0], xlon[jt,:], xlat[jt,:], pvIDs=[], cfig=cfig, cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1. )
            print('')
