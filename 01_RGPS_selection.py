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

#import climporn as cp
from climporn import chck4f, epoch2clock, clock2epoch

l_drop_doublons = True
rd_tol = 5. # tolerance distance in km to conclude it's the same buoy

# Time range of interest:
cdt1 = 'YYYY-01-01_00:00:00'
#
cdtI = 'YYYY-01-31_23:00:00' ; # "intermediate date": we drop buoys that do not go beyond this date (series would be too short)
#
cdt2 = 'YYYY-01-31_00:00:00'

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!

dt_buoy = 3*24*3600 ; # the expected mean time step of the input data, ~ 3 days [s]
dt_scan =    3*3600 ; # time increment while searching valid time intervals
dt_tolr =    1*3600 ; # time interval aka tolerance `+-dt_tolr` to consider two byoys are synchronized (Bouchat et al. 2021) [s]

idebug = 1 ; 

idbg_subsmpl = 0 ; vIDdbg = [ 18 , 62 , 98, 1200, 2800, 3000, 8800, 10290, 16805 ] ; l_show_IDs_fig = True ; # Debug: pick a tiny



list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

interp_1d = 0 ; # Time interpolation to fixed time axis: 0 => linear / 1 => akima



def plot_interp_series( iID, cname, vTs, vTt, vFs, vFt ):
    #
    # For debugging:
    # Overlay of original source F(t) series vFs(vTs)
    # and interpolated version vFt on vTt axis
    #
    import matplotlib.dates as mdates
    #
    cfig = 'debug_interp_'+cname+'_ID'+'%6.6i'%(iID)+'.'+fig_type
    fig = plt.figure(num=1, figsize=(12,5), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes([0.05, 0.13, 0.9, 0.8])
    #plt.axis([ min(vTt)-rdt, max(vTt)+rdt, min(xlat[:,ic])-0.01, max(xlat[:,ic])+0.01])
    #
    func = nmp.vectorize(dt.utcfromtimestamp)
    pl1 = plt.plot( mdates.date2num(func(vTs)), vFs, 'o-', color='#041a4d', linewidth=4., ms=10.)
    pl2 = plt.plot( mdates.date2num(func(vTt)), vFt, '*-', color='r'      , linewidth=1., ms=3)
    date_fmt = '%Y/%m/%d'
    date_formatter = mdates.DateFormatter(date_fmt)
    ax.xaxis.set_major_formatter(date_formatter)
    fig.autofmt_xdate()
    #
    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)






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

    rdt1 = clock2epoch(cdt1)
    rdtI = clock2epoch(cdtI)
    rdt2 = clock2epoch(cdt2)
    print( "   ===> in epoch time: ", rdt1, "to", rdt2 )
    print( "       ====> double check: ", epoch2clock(rdt1), "to",  epoch2clock(rdt2))


    # Build time axis at the buoys' `dt` for reference:
    rdt = dt_scan # time increment: 1d [s]

    Nts = int(round((rdt2 - rdt1) / rdt)) + 1
    print("\n *** New fixed time axis to use to scan data:")
    print( "   ===> Nts = "+str(Nts)+" days!")
    vTscan = nmp.zeros((Nts,3)) ; cTref = nmp.zeros(Nts, dtype='U19')
    vTscan[0,0] = rdt1 ; cTref[0] = epoch2clock(rdt1)
    for jt in range(1,Nts):
        vTscan[jt,0] = vTscan[jt-1,0] + rdt
        cTref[jt] = epoch2clock(vTscan[jt,0])
    # bounds:
    vTscan[:,1] = vTscan[:,0] - dt_tolr
    vTscan[:,2] = vTscan[:,0] + dt_tolr

        
    if idebug>0:
        for jt in range(Nts):
            print("   --- jt: "+str(jt)+" => ",vTscan[jt,0]," => ",cTref[jt])
            print("          => bounds: "+epoch2clock(vTscan[jt,1])+" - "+epoch2clock(vTscan[jt,2])+"\n")
            

    # Opening and inspecting input file
    chck4f(cf_in)

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

    vtime0 = nmp.zeros(Np0, dtype=int)
    vtime0 = id_in.variables['time'][:]

    vIDrgps0    = nmp.zeros(Np0, dtype=int)
    vIDrgps0[:] = id_in.variables['index'][:]

    id_in.close()

    
    vlon0[:] = nmp.mod(vlon0, 360.) ; # Longitudes in the [0:360] frame...


    # Masking all point that are before and beyond our period of interest:
    vmask = nmp.zeros(Np0, dtype=int) + 1
    vmask[nmp.where(vtime0 < rdt1-dt_tolr)] = 0
    vmask[nmp.where(vtime0 > rdt2-dt_tolr)] = 0

    vIDrgps0 =  nmp.ma.masked_where( vmask==0, vIDrgps0 )
    vtime0   =  nmp.ma.masked_where( vmask==0, vtime0   )
    vlon0    =  nmp.ma.masked_where( vmask==0, vlon0   )
    vlat0    =  nmp.ma.masked_where( vmask==0, vlat0   )
    

    
    if idbg_subsmpl>0:
        vIDs = nmp.array(vIDdbg, dtype=int)
    else:
        vIDs = nmp.sort( nmp.unique( vIDrgps0 ) )
        
    Nb   = len(vIDs)
    print("\n *** There are "+str(Nb)+" buoys to follow...")
    #exit(0)

    #
    Ns_max = 2000 # Max number of Stream

    XX = nmp.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the IDs used for a given stream...

    xstreams = []
    rt_prev_stream = 1e12
    istream = -1
    for jt in range(Nts):
        #
        rt = vTscan[jt,0]
        print("\n *** Selection of buoys that exist at "+cTref[jt]+" +-"+str(int(dt_tolr/3600))+"h!")
        #
        idx_ok, = nmp.where( nmp.abs( vtime0[:] - rt ) <= dt_tolr )
        Np      = len(idx_ok)
        #
        if idebug>1: print("    => "+str(Np)+" buoys satisfy this!")
        #
        if Np > 500 and rt < rt_prev_stream+dt_buoy -dt_tolr:
            #
            rt_prev_stream = rt
            #
            # That's a new stream
            istream = istream+1
            #
            if idebug>0: print("    => this date will the be the start of stream #"+str(istream)+", with "+str(Np)+" buoys!")
            
            
            # Which are the buoys:
            vidsT = vIDrgps0[idx_ok]
            nbb   = len(vidsT)

            ib = -1
            for jid in vidsT:
                ib = ib + 1
                idx_id, = nmp.where( vIDrgps0 == jid)
                vt = vtime0[idx_id]
                #
                if idebug>1:
                    print("    => ID="+str(jid)+": current and following times:")                    
                    for rt in vt:
                        print(epoch2clock(rt))
                #
                ntt = len(vt)
                if ntt > 1:
                    XX[istream,ib] = jid ; # keeps memory of buoys that have been used!
                    xstreams.append({
                        'istream': istream,
                        'id':      jid,
                        'Nrec':    ntt,
                        't0':      vt[0],
                        'tN':      vt[ntt-1],
                        'ct0':     epoch2clock(vt[0]),
                        'ctN':     epoch2clock(vt[ntt-1]),
                        })
                    #                    
        else:
            print("    => no buoys in this time range!")

        if istream >= 2: break ; #DEBUG!

    for js in range(2):
        print('\n *** stream #'+str(js)+' => IDs =')
        for jii in XX[js,:]:
            if jii>0: print(jii)



    exit(0)
    
    vt_earlst = nmp.zeros(Nb) ; ct_earlst = nmp.zeros(Nb, dtype='U20')
    vt_latest = nmp.zeros(Nb) ; ct_latest = nmp.zeros(Nb, dtype='U20')

    # For each buoy, getting date of earliest and latest snapshot
    print("\n *** Scanning to study precise date range of this file...")

    for jb in range(Nb):
        if jb%1000 == 0: print("    .... "+str(jb)+" / "+str(Nb)+"...")

        jid = vIDs[jb]

        idx = nmp.where( vIDrgps0==jid )
        vt  = vtime0[idx]
        t1  = nmp.min(vt) ; ct1 = epoch2clock(t1)
        t2  = nmp.max(vt) ; ct2 = epoch2clock(t2)

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
        print("  --- absolute earliest date is: "+epoch2clock(rt1))
        print("  --- absolute   latest date is: "+epoch2clock(rt2))
        print("  --- latest   of the earliests dates is: "+epoch2clock(rt1l))
        print("  --- earliest of the   latest  dates is: "+epoch2clock(rt2l),"\n")


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

        idx  = nmp.where( vIDrgps0==jid )
        vt   = vtime0[idx] ; # that's the time axis of this particular buoy [epoch time]
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
            if idebug>1: print("    * buoy with ID: "+str(jid)+" does not make it to "+cdt2, "\n      => dies at "+epoch2clock(rtend))
            for Ntz in range(Nt):
                rt = vTscan[Ntz]
                if rt > rtend: break
            if idebug>1: print("      ==> last index of vTscan to use is: "+str(Ntz-1)+" => "+epoch2clock(vTscan[Ntz-1]))

        #f = interpolate.interp1d(vt, vlat) ; vintrp = f(vTscan)
        #NUMPY: vintrp = nmp.interp(vTscan, vt, vlat) ; #, left=None, right=None, period=None)

        if   interp_1d==0:
            f = interpolate.interp1d(           vt, vlat0[idx])
        elif interp_1d==1:
            f = interpolate.Akima1DInterpolator(vt, vlat0[idx])
        xlat[0:Ntz,ic] = f(vTscan[0:Ntz])

        if   interp_1d==0:
            f = interpolate.interp1d(           vt, vlon0[idx])
        elif interp_1d==1:
            f = interpolate.Akima1DInterpolator(vt, vlon0[idx])
        xlon[0:Ntz,ic] = f(vTscan[0:Ntz])

        if l_shorten: xmsk[Ntz:Nt,ic] = 0

        if idebug>1 and Nb<=20:
            # Visual control of interpolation:
            plot_interp_series( jid, 'lat', vt, vTscan[0:Ntz], vlat0[idx], xlat[0:Ntz,ic] )
            plot_interp_series( jid, 'lon', vt, vTscan[0:Ntz], vlon0[idx], xlon[0:Ntz,ic] )


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
    from cartopy.crs import NorthPolarStereo
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
        v_time[jt]   = vTscan[jt]
        v_latb[jt,:] =  xlat[jt,:]
        v_lonb[jt,:] =  xlon[jt,:]
        v_y[jt,:]    =  yp[jt,:]
        v_x[jt,:]    =  xp[jt,:]

    f_out.About  = 'RGPS sea-ice drift data (Kwok, 1998)'
    f_out.Author = 'Generated with `'+path.basename(argv[0])+'` (L. Brodeau, 2022)'
    f_out.close()
    print("      ===> "+cf_out+" saved!")
