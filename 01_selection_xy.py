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
#

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np

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
#cdt2 = 'YYYY-01-31_00:00:00'

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!

dt_buoy = 3*24*3600 ; # the expected nominal time step of the input data, ~ 3 days [s]
dt_scan =    6*3600 ; # time increment while scanning for valid time intervals
#dt_scan =    12*3600 ; # time increment while scanning for valid time intervals
#dt_scan =    24*3600 ; # time increment while scanning for valid time intervals
#dt_scan =    36*3600 ; # time increment while scanning for valid time intervals
dt_tolr = dt_scan/2. ; # time interval aka tolerance `+-dt_tolr` to consider two byoys are synchronized (Bouchat et al. 2021) [s]

Ns_max = 200  # Max number of Streams, guess!!! #fixme...

Nb_min_stream = 200 ; # minimum number of buoys for considering a stream a stream!

Nb_min_cnsctv = 2   ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)

DistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]

idebug = 2 ;


list_expected_var = [ 'index', 'x', 'y', 'lon', 'lat', 'q_flag', 'time' ]

interp_1d = 0 ; # Time interpolation to fixed time axis: 0 => linear / 1 => akima



if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    if not path.exists('./npz'): mkdir('./npz')

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
    NTscan = int(round((rdt2 - rdt1) / dt_scan)) + 1
    print("\n *** New fixed time axis to use to scan data:")
    print( "   ===> NTscan = "+str(NTscan)+" days!")
    vTscan = np.zeros((NTscan,3), dtype=int  ) ; # `*,0` => precise time | `*,1` => bound below | `*,2` => bound above
    cTscan = np.zeros( NTscan   , dtype='U19')
    vTscan[0,0] =             rdt1
    cTscan[0]   = epoch2clock(rdt1)
    for jt in range(1,NTscan):
        tt = vTscan[jt-1,0] + dt_scan
        vTscan[jt,0] = tt
        cTscan[jt]   = epoch2clock(tt)
    # bounds:
    vTscan[:,1] = vTscan[:,0] - dt_tolr
    vTscan[:,2] = vTscan[:,0] + dt_tolr


    if idebug>0:
        for jt in range(NTscan):
            print("   --- jt: "+str(jt)+" => ",vTscan[jt,0]," => ",cTscan[jt])
            print("          => bounds: "+epoch2clock(vTscan[jt,1])+" - "+epoch2clock(vTscan[jt,2])+"\n")
    #exit(0)

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
        vtime0 = np.zeros(Np0, dtype=int)
        vtime0 = id_in.variables['time'][:]

        # Coordinates:
        vx0    = id_in.variables['x'][:]
        vy0    = id_in.variables['y'][:]
        vlon0  = id_in.variables['lon'][:]
        vlat0  = id_in.variables['lat'][:]

        # Buoys' IDs:
        vIDrgps0    = np.zeros(Np0, dtype=int)
        vIDrgps0[:] = id_in.variables['index'][:]


    vlon0[:] = np.mod(vlon0, 360.) ; # Longitudes in the [0:360] frame...


    # Masking all point that are before and beyond our period of interest:
    vmsk_time = np.zeros(Np0, dtype=int) + 1
    vmsk_time[np.where(vtime0 < rdt1-dt_tolr)] = 0
    vmsk_time[np.where(vtime0 > rdt2-dt_tolr)] = 0

    vIDrgps0 =  np.ma.masked_where( vmsk_time==0, vIDrgps0 )
    vtime0   =  np.ma.masked_where( vmsk_time==0, vtime0   )
    vx0      =  np.ma.masked_where( vmsk_time==0, vx0   )
    vy0      =  np.ma.masked_where( vmsk_time==0, vy0   )    
    vlon0    =  np.ma.masked_where( vmsk_time==0, vlon0   )
    vlat0    =  np.ma.masked_where( vmsk_time==0, vlat0   )

    # Remaining buoys (IDs)
    vIDs = np.sort( np.unique( vIDrgps0 ) )
    Nb   = len(vIDs)
    print("\n *** There are "+str(Nb)+" buoys to follow...")


    # Vectors along streams:
    VNB = [] ;  # Number of valid buoys in each stream

    # In the following, both Ns_max & Nb are excessive upper bound values #fixme
    XIDs = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the IDs used for a given stream...
    XNrc = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the number of records
    XTi  = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the actual initial time for the buoy
    XTe  = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # bad max size!! Stores the actual final time  "    "
    Xmsk = np.zeros((Ns_max, Nb), dtype=int)

    ID_in_use = []  ; # keeps memory of buoys that are already been considered!

    rT_prev_stream = 1e12
    istream        = -1
    for jt in range(NTscan):
        #
        rT = vTscan[jt,0] ; # current scan time
        #
        print("\n *** Selection of buoys that exist at "+cTscan[jt]+" +-"+str(int(dt_tolr/3600))+"h!")
        idx_ok, = np.where( np.abs( vtime0[:] - rT ) <= dt_tolr )
        Nok      = len(idx_ok)
        if idebug>1: print("    => "+str(Nok)+" buoys satisfy this!")

        Nbuoys_stream = 0
        if Nok > Nb_min_stream:
            # That's a new stream
            istream   = istream+1

            if idebug>0: print("    => this date will be the start of stream #"+str(istream)+", with "+str(Nok)+" buoys!")

            rT_prev_stream = rT
            
            vidsT = vIDrgps0[idx_ok] ; # IDs of the buoys that satisfy this

            ib = -1 ; # buoy counter...
            for jid in vidsT:
                #
                if not jid in ID_in_use:
                    #
                    ib = ib + 1
                    idx_id, = np.where( vIDrgps0 == jid)
                    #
                    vt1b  = vtime0[idx_id] ; # time records for this particular buoy
                    nbRec = len(vt1b)
                    #
                    if idebug>2:
                        print("    => ID="+str(jid)+": current and following times:")
                        for tt in vt1b: print(epoch2clock(tt))

                    # * Analysis of the time records for this particular buoy...
                    #    => Must cut off the series of the buoys as soon as its dt is too far
                    #       from the nominal time step:                    
                    #   => based on the initial time for this particular buoy:
                    #      - construct the ideal expected `vt1b` (`vt1b_ideal`) based on nominal dt
                    #      - dezing tout ce qui s'eloigne trop de ce vt1b_ideal !
                    vt1b_ideal = np.array( [ vt1b[0]+float(i)*float(dt_buoy) for i in range(nbRec) ], dtype=float )
                    print('LILO: vt1b =', vt1b)
                    print('LILO: vt1b_ideal =', vt1b_ideal)
                    lfu = np.any( np.abs(vt1b-vt1b_ideal)>dt_tolr)
                    print('LILO: fuckup? =>',lfu)
                    print('LILO')
                    #exit(0)
                    nbRecCut=nbRec
                    for jt in range(1,nbRec):
                        dt = vt1b[jt] - vt1b[jt-1]
                        print('dt = ',int(dt/3600))
                        terr = abs(dt - dt_buoy)
                        if terr > 0.05*dt_buoy:
                            print('LOLO CUTOFF!!! (buoy #'+str(jid)+') terr =',terr,'at jt = ',jt,'(nbRec=',nbRec,')')
                            nbRecCut = jt-1
                            break
                    print('   ===> nbRec, nbRecCut =',nbRec,nbRecCut,'\n')

                    #lilo!
                    #
                    # Initial position for the buoy:
                    it1, it2 = idx_id[0], idx_id[nbRec-1]
                    rd_ini = lbr.Dist2Coast( vlon0[it1], vlat0[it1], vlon_dist, vlat_dist, xdist )
                    rd_fin = lbr.Dist2Coast( vlon0[it2], vlat0[it2], vlon_dist, vlat_dist, xdist )
                    #print('\nLOLO: ==> initial distance to land =', rd_ini, 'km') ; exit(0)

                    # We want at least `Nb_min_cnsctv` consecutive records for the buoy
                    # and we want it to be located at least `DistFromLand` km off the coast
                    if nbRec >= Nb_min_cnsctv and rd_ini > DistFromLand and rd_fin > DistFromLand:

                        ID_in_use.append(jid)
                        Nbuoys_stream = Nbuoys_stream + 1 ; # this is another valid buoy for this stream
                        Xmsk[istream,ib] = 1      ; # valid point
                        XIDs[istream,ib] = jid    ; # keeps memory of buoys that have been used!
                        XNrc[istream,ib] = nbRec    ; # keeps memory of stream
                        #XTi[istream,ib] = vt1b[0]
                        #XTe[istream,ib] = vt1b[nbRec-1]
                        #
                    ### if nbRec >= Nb_min_cnsctv
                ### if not jid in ID_in_use
            ### for jid in vidsT

            if Nbuoys_stream>1:
                print("  * for now, we retained "+str(Nbuoys_stream)+" buoys in this stream...")
                VNB.append(Nbuoys_stream)
            else:
                print("  * well, none of the buoys of this stream had a at least "+str(Nb_min_cnsctv)+" following records!")
                istream = istream - 1
                if idebug>0: print("    => this was not a stream! So back to stream #"+str(istream)+" !!!")

    ### with Dataset(cf_in) as id_in:

    
    # Masking arrays:
    XIDs = np.ma.masked_where( Xmsk==0, XIDs )
    XNrc = np.ma.masked_where( Xmsk==0, XNrc )

    Nstreams = istream+1
    if len(VNB) != Nstreams: print('ERROR: number of streams?'); exit(0)
    print('\n *** Number of identified streams: '+str(Nstreams))


    for js in range(Nstreams):
        vids = np.ma.MaskedArray.compressed( XIDs[js,:] ) ; # valid IDs for current stream: shrinked, getting rid of masked points
        NvB  = VNB[js]
        if NvB != len(vids): print('ERROR Z1!'); exit(0)
        #
        print('\n\n *** Having a look at stream #'+str(js)+' !')
        print('     ===> has '+str(VNB[js])+' valid buoys!')
        if idebug>2:
            print('        => with following IDs:')
            for jii in vids: print(jii,' ', end="")
            print('')

        # Some scanning
        # 1- longest record of all the buoys
        Nt_max =  0
        jb_max = -1
        for jb in range(NvB):
            jid  = vids[jb]
            idx_id, = np.where( vIDrgps0 == jid) ; # => there can be only 2 (consecutive) points !!! See above!!!
            nr = len(idx_id)
            if nr>Nt_max:
                Nt_max = nr
                jb_max = jb
        Nt_lngst = np.max(XNrc[js,:]) ; # Max number of record from the buoy that has the most
        if Nt_lngst != Nt_max: print('ERROR: YY!'); exit(0)


        xx   = np.zeros((Nt_lngst,NvB))
        xy   = np.zeros((Nt_lngst,NvB))
        xlon = np.zeros((Nt_lngst,NvB))
        xlat = np.zeros((Nt_lngst,NvB))
        xtim = np.zeros((Nt_lngst,NvB) , dtype=int) -999 ; # the exact time for each buoy!
        xmsk = np.zeros((Nt_lngst,NvB) , dtype=int) + 1  ; # the mask for exluding buoys that stick out in time...


        
        # Show them on a map:
        for jb in range(NvB):
            jid  = vids[jb]
            idx_id, = np.where( vIDrgps0 == jid) ; # => there can be only 2 (consecutive) points !!! See above!!!
            nr = len(idx_id)
            #
            xx[:nr,jb]   =    vx0[idx_id]
            xy[:nr,jb]   =    vy0[idx_id]
            xlon[:nr,jb] =  vlon0[idx_id]
            xlat[:nr,jb] =  vlat0[idx_id]
            xtim[:nr,jb] = vtime0[idx_id]

        xmsk[np.where(xtim<=0)] = 0 ; # where time is zero, means that buoys does not exist anymore...
        xtim = np.ma.masked_where( xmsk==0, xtim ) ; # otherwize the `mean` in next line would use zeros!!!!
        vtim = np.mean(xtim, axis=1)

                
        # Nearest interpolation of vtim on the VTscan calendar !!!
        VT = vtim.copy()
        i=0
        for it in vtim:
            idx = np.argmin( np.abs(vTscan[:,0]-it) )
            if idebug>1: print('    it =',it,' => ',epoch2clock(it),' => nearest of VTscan =',epoch2clock(vTscan[idx,0]))
            VT[i] = vTscan[idx,0]
            i=i+1            
        del vtim
        
        # Are there buoys which are too far from this reference time ? lilo
        for jb in range(NvB):
            idt = np.abs(xtim[:,jb] - VT[:])
            if np.any(idt>dt_tolr):
                print('WOW, buoy #'+str(vids[jb])+' is more than '+str(int(dt_tolr/3600))
                      +' hours away from reference time...')
                (idw,) = np.where(idt>dt_tolr)
                print(idw)
                print(' ==> supressing '+str(len(idw))+' values!')
                xmsk[idw,jb] = 0
        #lilo

                        
        del xtim
        exit(0)



        
        # Saving 1 file per stream and per date:
        for jt in range(Nt_lngst):
            ct = split('_',epoch2clock(VT[jt]))[0] ; # just the day !
            cf_out = './npz/SELECTION_buoys_RGPS_stream'+'%3.3i'%(js)+'_'+ct+'.npz'
            np.savez_compressed( cf_out, itime=VT[jt], vx=xx[jt,:], vy=xy[jt,:], vlon=xlon[jt,:], vlat=xlat[jt,:], vids=vids[:] )

        if idebug>0:
            kf = lbr.ShowBuoysMap_Trec( VT, xlon, xlat, pvIDs=[], cnmfig='SELECTION_buoys_RGPS_stream'+'%3.3i'%(js), clock_res='d' )



