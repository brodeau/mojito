#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import argv, exit
from os import path
import numpy as nmp

from re import split

from netCDF4 import Dataset

from scipy import interpolate

from cartopy.crs import NorthPolarStereo

import climporn as cp
from   lbrgps import plot_interp_series

l_drop_doublons = True
rd_tol = 5. # tolerance distance in km to conclude it's the same buoy

# Time range of interest:
cdt1 = 'YYYY-01-01_00:00:00'
#
cdtI = 'YYYY-01-15_23:00:00' ; # "intermediate date": we drop buoys that do not go beyond this date (series would be too short)
#
cdt2 = 'YYYY-01-31_00:00:00'

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!


#idebug = 2 ; vIDdbg = [ 18 , 62 , 98, 1200, 2800, 3000, 8800, 10290, 16805 ] ; l_show_IDs_fig = True ; # Debug: pick a tiny selection of IDs to follow...
#idebug = 1 ; vIDdbg = range(0,36000,100)
idebug = 0

list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

interp_1d = 0 ; # Time interpolation to fixed time axis: 0 => linear / 1 => akima

#vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
#vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]





if __name__ == '__main__':

    narg = len(argv)
    if not narg in [3]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YEAR>')
        exit(0)
    cf_in = argv[1]
    cyear = argv[2]


    cdt1 = str.replace(cdt1, 'YYYY',cyear)
    cdtI = str.replace(cdtI, 'YYYY',cyear)
    cdt2 = str.replace(cdt2, 'YYYY',cyear)
    
    cf_out = 'RGPS_ice_drift_'+split('_', cdt1)[0]+'_'+split('_', cdt2)[0]+'_lb.nc' ;# netCDF file to generate


    print("\n *** Date range to restrain data to:")
    print(" ==> "+cdt1+" to "+cdt2 )
    print("      with buoys that do not pass "+cdtI+" canceled!")

    rdt1 = cp.clock2epoch(cdt1)
    rdtI = cp.clock2epoch(cdtI)
    rdt2 = cp.clock2epoch(cdt2)
    print( "   ===> in epoch time: ", rdt1, "to", rdt2 )
    print( "       ====> double check: ", cp.epoch2clock(rdt1), "to",  cp.epoch2clock(rdt2))


    # Build time axis we are going to interpolate onto:
    rdt = 3600.*24. # time increment: 1d [s]

    Nt = int(round((rdt2 - rdt1) / rdt)) + 1
    print("\n *** New fixed time axis to interpolate onto:")
    print( "   ===> Nt = "+str(Nt)+" days!")
    vTref = nmp.zeros(Nt) ; cTref = nmp.zeros(Nt, dtype='U19')
    vTref[0] = rdt1 ; cTref[0] = cp.epoch2clock(rdt1)
    for jt in range(1,Nt):
        vTref[jt] = vTref[jt-1] + rdt
        cTref[jt] = cp.epoch2clock(vTref[jt])
    #if idebug>0:
    for jt in range(Nt): print("   --- jt: "+str(jt)+" => ",vTref[jt]," => ",cTref[jt])


    # Opening and inspecting input file
    cp.chck4f(cf_in)

    id_in    = Dataset(cf_in)
    #
    list_dim = list( id_in.dimensions.keys() ) ; #print(' ==> dimensions: ', list_dim, '\n')
    if not 'points' in list_dim:
        print(' ERROR: no dimensions `points` found into input file!'); exit(0)
    #
    list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
    for cv in list_expected_var:
        if not cv in list_var:
            print(' ERROR: no variable `'+cv+'` found into input file!'); exit(0)
    #
    Np0 = id_in.dimensions['points'].size
    print(' *** Number of provided virtual buoys = ', Np0)


    vlon0  = id_in.variables['lon'][:]
    vlat0  = id_in.variables['lat'][:]

    rlat_min = nmp.min(vlat0)
    print('     ==> Southernmost latitude = ', rlat_min)

    rlat_max = nmp.max(vlat0)
    print('     ==> Northernmost latitude = ', rlat_max)

    ctunits = id_in.variables['time'].units
    if not ctunits == ctunits_expected:
        print(" ERROR: we expect '"+ctunits_expected+"' as units for time variables, yet we have: "+ctunits)

    vtime = nmp.zeros(Np0, dtype=int)
    vtime = id_in.variables['time'][:]

    vindex = nmp.zeros(Np0, dtype=int)
    vindex[:] = id_in.variables['index'][:]

    id_in.close()

    vlon0[:] = nmp.mod(vlon0, 360.) ; # Longitudes in the [0:360] frame...

    if idebug>0:
        vIDs = nmp.array(vIDdbg, dtype=int)
        cbla = "[DEBUG] IDs =>"
        for ii in vIDs: cbla = cbla+" "+str(ii)
    else:
        ## How many buoys (IDs) are there?
        id_min = nmp.min(vindex)
        id_max = nmp.max(vindex)
        vIDs = nmp.arange(id_min, id_max+1) ; # this should be all the buoys IDs
        cbla = "ID "+str(id_min)+" to ID "+str(id_max)+" !"

    Nb   = len(vIDs)
    print("\n *** There are "+str(Nb)+" buoys to follow: "+cbla)


    # First we are going to get rid of all buoys that are not fully spanning our defined period
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # To mask IDs we cancel:
    vmask = nmp.zeros(Nb, dtype=int)
    vmask[:] = 1

    vt_earlst = nmp.zeros(Nb) ; ct_earlst = nmp.zeros(Nb, dtype='U20')
    vt_latest = nmp.zeros(Nb) ; ct_latest = nmp.zeros(Nb, dtype='U20')

    # For each buoy, getting date of earliest and latest snapshot
    print("\n *** Scanning to study precise date range of this file...")

    for jb in range(Nb):
        if jb%1000 == 0: print("    .... "+str(jb)+" / "+str(Nb)+"...")

        jid = vIDs[jb]

        idx = nmp.where( vindex==jid )
        vt  = vtime[idx]
        t1  = nmp.min(vt) ; ct1 = cp.epoch2clock(t1)
        t2  = nmp.max(vt) ; ct2 = cp.epoch2clock(t2)

        # Get rid of IDs that are seeded after rdt1:
        if t1>rdt1:
            print("   --- excluding buoy #"+str(jb)+" with ID: ",jid," (pops up after "+cdt1+"!)")
            vmask[jb] = 0
        if t2<=rdtI:
            print("   --- excluding buoy #"+str(jb)+" with ID: ",jid," (vanishes before "+cdtI+"!)")
            vmask[jb] = 0

        if idebug>0 and Nb<=20: print("      * buoy #"+str(jb)+" (id="+str(jid)+": "+ct1+" ==> "+ct2)
        vt_earlst[jb] = t1 ; ct_earlst[jb] = ct1 ;
        vt_latest[jb] = t2 ; ct_latest[jb] = ct2 ;


    NbR = nmp.sum(vmask)
    print("\n *** Only "+str(NbR)+" buoys out of "+str(Nb)+" are remaining (time range)")

    # Masking all destroyed buoys:
    vIDs      =  nmp.ma.masked_where( vmask==0, vIDs      )
    vt_earlst =  nmp.ma.masked_where( vmask==0, vt_earlst )
    vt_latest =  nmp.ma.masked_where( vmask==0, vt_latest )

    if idebug>1:
        rt1  = nmp.min(vt_earlst)
        rt2  = nmp.max(vt_latest)
        rt1l = nmp.max(vt_earlst)
        rt2l = nmp.min(vt_latest)
        print("  --- absolute earliest date is: "+cp.epoch2clock(rt1))
        print("  --- absolute   latest date is: "+cp.epoch2clock(rt2))
        print("  --- latest   of the earliests dates is: "+cp.epoch2clock(rt1l))
        print("  --- earliest of the   latest  dates is: "+cp.epoch2clock(rt2l),"\n")


    # Dropping canceled (masked) buoys and updating number of valid buoys:
    vIDs = nmp.ma.MaskedArray.compressed(vIDs)
    Nb   = len(vIDs)
    print("\n *** UPDATE: there are now "+str(Nb)+" buoys to follow!")


    # Final arrays have 2 dimmensions Nb & Nt (buoy ID & time record)
    xmsk = nmp.zeros((Nt,Nb), dtype=int)
    xlon = nmp.zeros((Nt,Nb))
    xlat = nmp.zeros((Nt,Nb))

    xmsk[:,:] = 1

    ic = -1
    for jid in vIDs:
        ic = ic + 1

        #if vmask[ic] == 1:

        idx  = nmp.where( vindex==jid )
        vt   = vtime[idx] ; # that's the time axis of this particular buoy [epoch time]
        Np   = len(vt)    ; # Mind that Np varies from 1 buoy to another...

        for jp in range(1,Np):
            if vt[jp]<=vt[jp-1]:
                print("ERROR: time not increasing","\n   ==> jp =",jp, " ==> ", vt[jp-1], " -- ", vt[jp])
                exit(0)

        Ntz = Nt

        # Does this buoy survive beyond the fixed time range or not??? lilo
        rtend = nmp.max(vt)
        l_shorten = ( rtend < rdt2 )
        if l_shorten:
            if idebug>1: print("    * buoy with ID: "+str(jid)+" does not make it to "+cdt2, "\n      => dies at "+cp.epoch2clock(rtend))
            for Ntz in range(Nt):
                rt = vTref[Ntz]
                if rt > rtend: break
            if idebug>1: print("      ==> last index of vTref to use is: "+str(Ntz-1)+" => "+cp.epoch2clock(vTref[Ntz-1]))

        #f = interpolate.interp1d(vt, vlat) ; vintrp = f(vTref)
        #NUMPY: vintrp = nmp.interp(vTref, vt, vlat) ; #, left=None, right=None, period=None)

        if   interp_1d==0:
            f = interpolate.interp1d(           vt, vlat0[idx])
        elif interp_1d==1:
            f = interpolate.Akima1DInterpolator(vt, vlat0[idx])
        xlat[0:Ntz,ic] = f(vTref[0:Ntz])

        if   interp_1d==0:
            f = interpolate.interp1d(           vt, vlon0[idx])
        elif interp_1d==1:
            f = interpolate.Akima1DInterpolator(vt, vlon0[idx])
        xlon[0:Ntz,ic] = f(vTref[0:Ntz])

        if l_shorten: xmsk[Ntz:Nt,ic] = 0

        if idebug>1 and Nb<=20:
            # Visual control of interpolation:
            lbr.plot_interp_series( jid, 'lat', vt, vTref[0:Ntz], vlat0[idx], xlat[0:Ntz,ic] )
            lbr.plot_interp_series( jid, 'lon', vt, vTref[0:Ntz], vlon0[idx], xlon[0:Ntz,ic] )


    #
    xlat = nmp.ma.masked_where( xmsk==0, xlat )
    xlon = nmp.ma.masked_where( xmsk==0, xlon )

    
    # Removing doublons:
    if l_drop_doublons:
        from gonzag import Haversine

        vid_mask = nmp.zeros(Nb, dtype=int)
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
                if idebug>0: print("\n *** 'TOO CLOSE' scanning: Buoy #"+str(jb)+"=> ID ="+str(jid))
        
                rlat = vlat[jb]
                rlon = vlon[jb]
        
                if idebug>0: print("    ==> lat,lon = ",rlat,rlon)
                
                # Building array of distances (km) with all other buoy at this same time record:
                vdist = Haversine( rlat, rlon, vlat, vlon )
        
                # Need to mask itself (distance = 0!)
                if vdist[jb] != 0.:
                    print(" PROBLEM: distance with yourself should be 0!")
                vdist[jb] = 9999.
        
                if idebug>1:
                    for rd in vdist:
                        if rd<10.: print("   distance = ", rd)
        
                rdmin = nmp.min(vdist)
                if idebug>0: print("    ==> closest other buoy is at "+str(round(rdmin,2))+" km !")
                
                if rdmin < rd_tol:
                    # There is at least 1 buoy too close!
                    #  => populating the buoys closer than `rd_tol` km:
                    (vidx_close,) = nmp.where( vdist < rd_tol )
        
                    if idebug>0:
                        print("    ==> list of buoys closer than "+str(rd_tol)+" km:")
                        for jb in vidx_close:
                            print(vIDs[ii])
    
                    vid_mask[(vidx_close,)] = 0

        # Alright, all IDs spotted via vid_mask should be cancelled in the 2D fields, at all time steps
        (idmsk,) = nmp.where(vid_mask==0)
        xmsk[:,idmsk] = 0
        vIDs = nmp.ma.masked_where( vid_mask==0, vIDs )
        xlat = nmp.ma.masked_where( xmsk==0, xlat )
        xlon = nmp.ma.masked_where( xmsk==0, xlon )
        
        Nok = nmp.sum(xmsk[0,:])
        print("\n *** UPDATE: after 'too close' canceling, there are "+str(Nok)+" buoys left at t0!")

        # Need to compress arrays, didn't find an elegant way to do it so here is the following abomination:
        xmsk2 = nmp.zeros((Nt,Nok), dtype=int)
        xlon2 = nmp.zeros((Nt,Nok))
        xlat2 = nmp.zeros((Nt,Nok))
        ic = 0
        for jb in range(Nb):
            if xmsk[0,jb] == 1:
                xmsk2[:,ic] = xmsk[:,jb]
                xlon2[:,ic] = xlon[:,jb]
                xlat2[:,ic] = xlat[:,jb]
                ic = ic+1
        vIDs = nmp.ma.MaskedArray.compressed(vIDs)
        Nb = Nok
        #print("Shape of vIDs, xlat2, xlon2, xmsk2 =",vIDs.shape,xlat2.shape,xlon2.shape,xmsk2.shape)
        del xmsk, xlat, xlon
        xmsk = xmsk2
        xlat = xlat2
        xlon = xlon2
        del xmsk2, xlat2, xlon2
        #print("Shape of vIDs, xlat, xlon, xmsk =",vIDs.shape,xlat.shape,xlon.shape,xmsk.shape)
        #exit(0)
        
    #endif l_drop_doublons

    # convert to meters in neXtSIM projection    
    srs_nextsim = NorthPolarStereo(central_longitude=-45, true_scale_latitude=60)
    srs_rgps    = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
    xpos,ypos,_ = srs_nextsim.transform_points(srs_rgps, xlon, xlat).T * 1000.

    xp = xpos.T
    yp = ypos.T
    del xpos,ypos

    yp   = nmp.ma.masked_where( xmsk==0, yp )
    xp   = nmp.ma.masked_where( xmsk==0, xp )
    xlat = nmp.ma.masked_where( xmsk==0, xlat )
    xlon = nmp.ma.masked_where( xmsk==0, xlon )

    
    
    
    # GENERATION OF COMPREHENSIVE NETCDF FILE
    #########################################


    print("\n *** About to generate file: "+cf_out+" ...")

    f_out = Dataset(cf_out, 'w', format='NETCDF4')

    cd_time = 'time'
    cd_buoy = 'id_buoy'

    #Dimensions:
    f_out.createDimension(cd_time, None)
    f_out.createDimension(cd_buoy, Nb  )

    v_time  = f_out.createVariable(cd_time, 'i4',(cd_time,))
    v_id    = f_out.createVariable(cd_buoy, 'i4',(cd_buoy,))
    v_latb  = f_out.createVariable('latitude' , 'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)
    v_lonb  = f_out.createVariable('longitude', 'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)
    v_y     = f_out.createVariable('y_pos' ,    'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)
    v_x     = f_out.createVariable('x_pos',     'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)

    v_time.units = ctunits
    v_id.units   = "ID of buoy"
    v_latb.units = "degrees north"
    v_lonb.units = "degrees south"
    v_y.units    = "m"
    v_x.units    = "m"

    v_id[:] = vIDs[:]

    for jt in range(Nt):
        v_time[jt]   = vTref[jt]
        v_latb[jt,:] =  xlat[jt,:]
        v_lonb[jt,:] =  xlon[jt,:]
        v_y[jt,:]    =  yp[jt,:]
        v_x[jt,:]    =  xp[jt,:]

    f_out.About  = 'RGPS sea-ice drift data (Kwok, 1998)'
    f_out.Author = 'Generated with `'+path.basename(argv[0])+'` (L. Brodeau, 2022)'
    f_out.close()
    print("      ===> "+cf_out+" saved!")
