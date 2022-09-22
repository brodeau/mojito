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

from climporn import chck4f, epoch2clock, clock2epoch

import lbrgps as lbr


#l_drop_doublons = True
#rd_tol = 5. # tolerance distance in km to conclude it's the same buoy

# Time range of interest:
cdt1 = 'YYYY-01-01_00:00:00'
cdt2 = 'YYYY-01-10_00:00:00'

fdist2coast_nc = '/MEDIA/data/data/dist2coast/dist2coast_4deg_North.nc'

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!

dt_buoy = 3*24*3600 ; # the expected nominal time step of the input data, ~ 3 days [s]
dt_scan =    3*3600 ; # time increment while scanning for valid time intervals
dt_tolr =    1*3600 ; # time interval aka tolerance `+-dt_tolr` to consider two byoys are synchronized (Bouchat et al. 2021) [s]

Ns_max = 200  # Max number of Streams, guess!!! #fixme...

Nb_min_stream = 100 ; # minimum number of buoys for considering a stream a stream!

Nb_min_cnsctv = 2   ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)

DistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]

idebug = 2 ;

#idbg_subsmpl = 0 ; vIDdbg = [ 18 , 62 , 98, 1200, 2800, 3000, 8800, 10290, 16805 ] ; l_show_IDs_fig = True ; # Debug: pick a tiny



list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

interp_1d = 0 ; # Time interpolation to fixed time axis: 0 => linear / 1 => akima



if __name__ == '__main__':

    #print('LOLO: rd =', rd, 'km !')
    #exit(0)

    narg = len(argv)
    if not narg in [3]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YEAR>')
        exit(0)
    cf_in = argv[1]
    cyear = argv[2]


    cdt1 = str.replace(cdt1, 'YYYY',cyear)
    cdt2 = str.replace(cdt2, 'YYYY',cyear)

    cf_out = 'RGPS_ice_drift_'+split('_', cdt1)[0]+'_'+split('_', cdt2)[0]+'_lb.nc' ;# netCDF file to generate


    print("\n *** Date range to restrain data to:")
    print(" ==> "+cdt1+" to "+cdt2 )

    rdt1 = clock2epoch(cdt1)
    rdt2 = clock2epoch(cdt2)
    print( "   ===> in epoch time: ", rdt1, "to", rdt2 )
    print( "       ====> double check: ", epoch2clock(rdt1), "to",  epoch2clock(rdt2))


    # Load `distance to coast` data:
    vlon_dist, vlat_dist, xdist = lbr.LoadDist2CoastNC( fdist2coast_nc )

    # Build scan time axis willingly at relative high frequency (dt_scan << dt_buoy)
    Nts = int(round((rdt2 - rdt1) / dt_scan)) + 1
    print("\n *** New fixed time axis to use to scan data:")
    print( "   ===> Nts = "+str(Nts)+" days!")
    vTscan = nmp.zeros((Nts,3), dtype=int  ) ; # `*,0` => precise time | `*,1` => bound below | `*,2` => bound above
    cTscan = nmp.zeros( Nts   , dtype='U19')
    vTscan[0,0] =             rdt1
    cTscan[0]   = epoch2clock(rdt1)
    for jt in range(1,Nts):
        tt = vTscan[jt-1,0] + dt_scan
        vTscan[jt,0] = tt
        cTscan[jt]   = epoch2clock(tt)
    # bounds:
    vTscan[:,1] = vTscan[:,0] - dt_tolr
    vTscan[:,2] = vTscan[:,0] + dt_tolr


    if idebug>0:
        for jt in range(Nts):
            print("   --- jt: "+str(jt)+" => ",vTscan[jt,0]," => ",cTscan[jt])
            print("          => bounds: "+epoch2clock(vTscan[jt,1])+" - "+epoch2clock(vTscan[jt,2])+"\n")


    # Opening and inspecting the input file
    #######################################
    chck4f(cf_in)

    with Dataset(cf_in) as id_in:

        list_dim = list( id_in.dimensions.keys() ) ; #print(' ==> dimensions: ', list_dim, '\n')
        if not 'points' in list_dim:
            print(' ERROR: no dimensions `points` found into input file!'); exit(0)

        list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
        for cv in list_expected_var:
            if not cv in list_var:
                print(' ERROR: no variable `'+cv+'` found into input file!'); exit(0)

        Np0 = id_in.dimensions['points'].size
        print(' *** Number of provided virtual buoys = ', Np0)

        # Time records:
        ctunits = id_in.variables['time'].units
        if not ctunits == ctunits_expected:
            print(" ERROR: we expect '"+ctunits_expected+"' as units for the time record vector, yet we have: "+ctunits)
            exit(0)
        vtime0 = nmp.zeros(Np0, dtype=int)
        vtime0 = id_in.variables['time'][:]

        # Coordinates:
        vlon0  = id_in.variables['lon'][:]
        vlat0  = id_in.variables['lat'][:]
        #rlat_min = nmp.min(vlat0)
        #print('     ==> Southernmost latitude = ', rlat_min)
        #rlat_max = nmp.max(vlat0)
        #print('     ==> Northernmost latitude = ', rlat_max)

        # Buoys' IDs:
        vIDrgps0    = nmp.zeros(Np0, dtype=int)
        vIDrgps0[:] = id_in.variables['index'][:]


    vlon0[:] = nmp.mod(vlon0, 360.) ; # Longitudes in the [0:360] frame...


    # Masking all point that are before and beyond our period of interest:
    vmsk_time = nmp.zeros(Np0, dtype=int) + 1
    vmsk_time[nmp.where(vtime0 < rdt1-dt_tolr)] = 0
    vmsk_time[nmp.where(vtime0 > rdt2-dt_tolr)] = 0

    vIDrgps0 =  nmp.ma.masked_where( vmsk_time==0, vIDrgps0 )
    vtime0   =  nmp.ma.masked_where( vmsk_time==0, vtime0   )
    vlon0    =  nmp.ma.masked_where( vmsk_time==0, vlon0   )
    vlat0    =  nmp.ma.masked_where( vmsk_time==0, vlat0   )

    # Remaining buoys (IDs)
    vIDs = nmp.sort( nmp.unique( vIDrgps0 ) )
    Nb   = len(vIDs)
    print("\n *** There are "+str(Nb)+" buoys to follow...")


    # Vectors along streams:
    VNB = [] ;  # Number of valid buoys in each stream
    VTi = [] ;  # Nominal initial time we retain for each stream

    # In the following, both Ns_max & Nb are excessive upper bound values #fixme
    XIDs = nmp.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the IDs used for a given stream...
    XNrc = nmp.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the number of records
    XT1  = nmp.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the number of records
    XT2  = nmp.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the number of records
    Xmsk = nmp.zeros((Ns_max, Nb), dtype=int)


    ID_in_use = []  ; # keeps memory of buoys that are already been considered!
    xstreams  = []  ; # array of dictionnary with informations to store...

    rt_prev_stream = 1e12
    istream        = -1
    for jt in range(Nts):
        #
        rt = vTscan[jt,0] ; # current scan time
        #
        print("\n *** Selection of buoys that exist at "+cTscan[jt]+" +-"+str(int(dt_tolr/3600))+"h!")
        idx_ok, = nmp.where( nmp.abs( vtime0[:] - rt ) <= dt_tolr )
        Nok      = len(idx_ok)
        if idebug>1: print("    => "+str(Nok)+" buoys satisfy this!")
        #
        Nbuoys_stream = 0

        #if Nok > Nb_min_stream and rt < rt_prev_stream+dt_buoy -dt_tolr:
        if Nok > Nb_min_stream :
            #
            # That's a new stream
            istream   = istream+1
            #
            if idebug>0: print("    => this date will be the start of stream #"+str(istream)+", with "+str(Nok)+" buoys!")

            rt_prev_stream = rt

            # Which are the buoys:
            vidsT = vIDrgps0[idx_ok]

            ib = -1 ; # buoy counter...

            for jid in vidsT:
                #
                if not jid in ID_in_use:
                    #
                    ib = ib + 1
                    idx_id, = nmp.where( vIDrgps0 == jid)
                    #
                    vt  = vtime0[idx_id] ; # time records for this particular buoy
                    ntt = len(vt)
                    #
                    if idebug>2:
                        print("    => ID="+str(jid)+": current and following times:")
                        for rt in vt: print(epoch2clock(rt))
                    #
                    #
                    # Initial position for the buoy:
                    tlon = vlon0[idx_id[0]]
                    tlat = vlat0[idx_id[0]]
                    rd = lbr.Dist2Coast( tlon, tlat, vlon_dist, vlat_dist, xdist )
                    #print('\nLOLO: ==> distance to land =', rd, 'km')
                    #exit(0)

                    # We want at least `Nb_min_cnsctv` consecutive records for the buoy
                    # and we want it to be located at least `DistFromLand` km off the coast
                    if ntt >= Nb_min_cnsctv and rd > DistFromLand:

                        ID_in_use.append(jid)
                        Nbuoys_stream = Nbuoys_stream + 1 ; # this is another valid buoy for this stream
                        Xmsk[istream,ib] = 1      ; # valid point
                        XIDs[istream,ib] = jid    ; # keeps memory of buoys that have been used!
                        XNrc[istream,ib] = ntt    ; # keeps memory of stream
                        XT1[ istream,ib] = vt[0]  ; # keeps memory of
                        XT2[ istream,ib] = vt[ntt-1] ; # keeps memory of
                        #
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
                    ### if ntt >= Nb_min_cnsctv
                ### if not jid in ID_in_use
            ### for jid in vidsT

            if Nbuoys_stream>1:
                print("  * we retained "+str(Nbuoys_stream)+" buoys in this stream...")
                VNB.append(Nbuoys_stream)
                VTi.append(rt)
            else:
                print("  * well, none of the buoys of this stream had a at least "+str(Nb_min_cnsctv)+" following records!")
                istream = istream - 1
                if idebug>0: print("    => this was not a stream! So back to stream #"+str(istream)+" !!!")


    # Masking arrays:
    XIDs = nmp.ma.masked_where( Xmsk==0, XIDs )
    XIDs = nmp.ma.masked_where( Xmsk==0, XIDs )

    Nstreams = istream+1
    if len(VNB) != Nstreams: print('ERROR: number of streams?'); exit(0)
    print('\n *** Number of identified streams: '+str(Nstreams))


    if idebug>0:
        for js in range(Nstreams):
            vids = nmp.ma.MaskedArray.compressed( XIDs[js,:] ) ; # valid IDs for current stream: shrinked, getting rid of masked points
            Nvb  = VNB[js]
            if Nvb != len(vids): print('ERROR Z1!'); exit(0)
            #
            print('\n\n *** Having a look at stream #'+str(js)+' !')
            print('     ===> has '+str(VNB[js])+' valid buoys!')
            if idebug>1:
                print('        => with following IDs:')
                for jii in vids: print(jii,' ', end="")
                print('')
            #

            # Some scanning
            # 1- longest record of all the byous
            Nt_max =  0
            jb_max = -1
            for jb in range(Nvb):
                jid  = vids[jb]
                idx_id, = nmp.where( vIDrgps0 == jid) ; # => there can be only 2 (consecutive) points !!! See above!!!
                nr = len(idx_id)
                if nr>Nt_max:
                    Nt_max = nr
                    jb_max = jb
            Nt_lngst = nmp.max(XNrc[js,:]) ; # Max number of record from the buoy that has the most
            if Nt_lngst != Nt_max: print('ERROR: YY!'); exit(0)

            vt   = nmp.zeros( Nt_lngst , dtype=int)
            xlon = nmp.zeros((Nt_lngst,Nvb))
            xlat = nmp.zeros((Nt_lngst,Nvb))

            # Show them on a map:
            for jb in range(Nvb):
                jid  = vids[jb]
                idx_id, = nmp.where( vIDrgps0 == jid) ; # => there can be only 2 (consecutive) points !!! See above!!!
                #if idebug>1: print("ID = "+str(jid)+" => points =>",idx_id)
                nr = len(idx_id)
                if jb==jb_max: vt[:] = vtime0[idx_id]
                xlon[:nr,jb] = vlon0[idx_id]
                xlat[:nr,jb] = vlat0[idx_id]



            #kf = lbr.ShowBuoysMap( VTi[js], vlon[:], vlat[:], pvIDs=vids, cnmfig='SELECTION_buoys_RGPS' )
            kf = lbr.ShowBuoysMap_Trec( vt, xlon, xlat, pvIDs=[], cnmfig='SELECTION_buoys_RGPS_stream'+'%3.3i'%(js) )
        exit(0)


        #    # Show them on a map:
        #    vlon = []
        #    vlat = []
        #    for jb in range(Nvb):
        #        jid  = vids[jb]
        #        idx_id, = nmp.where( vIDrgps0 == jid) ; # => there can be only 2 (consecutive) points !!! See above!!!
        #        if idebug>2: print("ID = "+str(jid)+" => points =>",idx_id)
        #        ip = idx_id[0] ; # Only initial point!!!
        #        #print(ip); exit(0)
        #        vlon.append( vlon0[ip] )
        #        vlat.append( vlat0[ip] )
        #
        #    #kf = lbr.ShowBuoysMap( VTi[js], vlon[:], vlat[:], pvIDs=vids, cnmfig='SELECTION_buoys_RGPS' )
        #    kf = lbr.ShowBuoysMap( VTi[js], vlon[:], vlat[:], pvIDs=[], cnmfig='SELECTION_buoys_RGPS_stream'+'%3.3i'%(js) )





