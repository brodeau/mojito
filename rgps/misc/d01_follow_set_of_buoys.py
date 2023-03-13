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

'''

   Will follow a sub-set of buoys during their whole life

'''


from sys import argv, exit
from os import path, environ, mkdir
import numpy as np
from re import split
from netCDF4 import Dataset

import mojito as mjt

idebug = 2

frqFollow = 400 ; # subsampling for the buoys to follow!

cdt_pattern = 'YYYY-MM-DD_00:00:00' ; # pattern for dates

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

dt_buoy = 3*24*3600 ; # the expected nominal time step of the input data, ~ 3 days [s]
dt_scan =    6*3600 ; # time increment while scanning for valid time intervals
dt_tolr = dt_scan/2. ; # time interval aka tolerance `+-dt_tolr` to consider two byoys are synchronized (Bouchat et al. 2021) [s]

Ns_max = 200  # Max number of Streams, guess!!! #fixme...

Nb_min_stream = 200 ; # minimum number of buoys for considering a stream a stream!

Nb_min_cnsctv = 2   ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)

MinDistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]

Nb_min_buoys = 50 ; # minimum number of buoys necessary for considering the record of a stream a record!

list_expected_var = [ 'index', 'x', 'y', 'lon', 'lat', 'q_flag', 'time' ]

interp_1d = 0 ; # Time interpolation to fixed time axis: 0 => linear / 1 => akima

FillValue = -9999. 

if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    for cd in ["npz","figs"]:
        if not path.exists('./'+cd): mkdir('./'+cd)    
    if not path.exists('./figs/SELECTION'): mkdir('./figs/SELECTION')    
        
    narg = len(argv)
    if not narg in [5]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YEAR> <MMDD1> <MMDD2>')
        exit(0)
    cf_in = argv[1]
    cyear = argv[2]
    cmmd1 = argv[3]
    cmmd2 = argv[4]

    cmm1, cdd1 = cmmd1[0:2], cmmd1[2:4]
    cmm2, cdd2 = cmmd2[0:2], cmmd2[2:4]
    cdt1 = str.replace(cdt_pattern,'YYYY',cyear) ; cdt1 = str.replace(cdt1,'MM',cmm1) ; cdt1 = str.replace(cdt1,'DD',cdd1)
    cdt2 = str.replace(cdt_pattern,'YYYY',cyear) ; cdt2 = str.replace(cdt2,'MM',cmm2) ; cdt2 = str.replace(cdt2,'DD',cdd2)
    
    print("\n *** Date range to restrain data to:")
    print(" ==> "+cdt1+" to "+cdt2 )
    
    rdt1, rdt2 = mjt.clock2epoch(cdt1), mjt.clock2epoch(cdt2)
    print( "   ===> in epoch time: ", rdt1, "to", rdt2 )
    print( "       ====> double check: ", mjt.epoch2clock(rdt1), "to",  mjt.epoch2clock(rdt2))


    # Opening and inspecting the input file
    #######################################
    mjt.chck4f(cf_in)

    with Dataset(cf_in) as id_in:

        list_dim = list( id_in.dimensions.keys() ) ; #print(' ==> dimensions: ', list_dim, '\n')
        if not 'points' in list_dim:
            print(' ERROR: no dimensions `points` found into input file!'); exit(0)

        list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
        for cv in list_expected_var:
            if not cv in list_var:
                print(' ERROR: no variable `'+cv+'` found into input file!'); exit(0)

        Np0 = id_in.dimensions['points'].size
        print('\n *** Total number of points in the file = ', Np0)

        # Time records:
        ctunits = id_in.variables['time'].units
        if not ctunits == mjt.tunits_default:
            print(" ERROR: we expect '"+mjt.tunits_default+"' as units for the time record vector, yet we have: "+ctunits)
            exit(0)
        vtime0 = id_in.variables['time'][:]

        # Coordinates:
        #vx0    = id_in.variables['x'][:]
        #vy0    = id_in.variables['y'][:]
        vlon0  = id_in.variables['lon'][:]
        vlat0  = id_in.variables['lat'][:]

        # Buoys' IDs:
        vBIDs0    = np.zeros(Np0, dtype=int) - 1
        vBIDs0[:] = id_in.variables['index'][:]

    ### with Dataset(cf_in) as id_in
    
    vlon0[:] = np.mod(vlon0, 360.) ; # Longitudes in the [0:360] frame...

    # Masking all point that are before and beyond our period of interest:
    vmsk_time = np.zeros(Np0, dtype=int) + 1
    vmsk_time[np.where(vtime0 < rdt1-dt_tolr)] = 0
    vmsk_time[np.where(vtime0 > rdt2+dt_tolr)] = 0
    #
    (idx_masked,) = np.where( vmsk_time == 0 )

    if Np0-len(idx_masked) != np.sum(vmsk_time):
        print('ERROR: fuck up #1!')
        exit(0)
    
    print('\n *** Total number of points remaining after time-range-exclusion = ',Np0-len(idx_masked), '=', np.sum(vmsk_time))
    #
    #
    vBIDs0[idx_masked] = int(FillValue) ; vBIDs0 =  np.ma.masked_where( vmsk_time==0, vBIDs0 )
    vtime0[idx_masked] = FillValue      ; vtime0 =  np.ma.masked_where( vmsk_time==0, vtime0 )
    #vx0[idx_masked]    = FillValue      ; vx0    =  np.ma.masked_where( vmsk_time==0, vx0    )
    #vy0[idx_masked]    = FillValue      ; vy0    =  np.ma.masked_where( vmsk_time==0, vy0    )
    vlon0[idx_masked]  = FillValue      ; vlon0  =  np.ma.masked_where( vmsk_time==0, vlon0  )
    vlat0[idx_masked]  = FillValue      ; vlat0  =  np.ma.masked_where( vmsk_time==0, vlat0  )

    # Remaining buoys (IDs)
    (idx,) = np.where(vBIDs0.data >= 0)
    vIDs = np.sort( np.unique( vBIDs0[idx] ) ) ; # if not `[idx]` then `FillValue` is counted once!!!
    Nb   = len(vIDs)
    print("\n *** We found "+str(Nb)+" different buoys alive during specified period of time!")

    IDs2Follow = []
    ic=0
    for jid in vIDs:
        #
        if ic%frqFollow == 0: IDs2Follow.append(jid)
        #
        ic = ic+1
    ### for jid in vIDs
    
    IDs2Follow = np.array(IDs2Follow, dtype=int)

    nBf = len(IDs2Follow)
    print('  ==> with a sub-sambpling of '+str(frqFollow)+', we have now '+str(nBf)+' to follow...')


    # Maximum number of records in all the buoys?
    nRmax = 0
    id_longest_rec = 0
    for jid in IDs2Follow:        
        (idx_id,) = np.where( vBIDs0 == jid)
        nR = max( len(idx_id), nRmax )

        if nR >= nRmax:
            nRmax = nR
            id_longest_rec = jid
        #print(idx_id)
        #for ii in idx_id:
        #    print(' * time: ',mjt.epoch2clock( vtime0[ii]) )

    print('\n *** The buoy with most records is '+str(id_longest_rec)+', it has '+str(nRmax)+' of them.')

    xmsk = np.zeros((nBf,nRmax), dtype=int) + 1
    xlon = np.zeros((nBf,nRmax))
    xlat = np.zeros((nBf,nRmax))
    xtim = np.zeros((nBf,nRmax))
    valr = np.zeros( nBf,        dtype=int) ; # number of valid records of each buoy

    ib = 0
    for jid in IDs2Follow:        
        (idx_id,) = np.where( vBIDs0 == jid)
        nr = len(idx_id)
        valr[ib]      = nr
        xlon[ib,0:nr] = vlon0[idx_id]
        xlat[ib,0:nr] = vlat0[idx_id]
        xtim[ib,0:nr] = vtime0[idx_id]
        #        
        xmsk[ib,nr:] = 0 ; # Following records do not exist for this buoy...
        #
        # Is time for this buoys monotonically increasing along records?
        if not ( np.all(np.diff(xtim[ib,0:nr]) > 0.) ):
            print('PROBLEM: time is not monotonically increasing along records for buoy with ID:'+str(jid)+' !!!')
            exit(0)
        #
        ib = ib+1

    xlon = np.ma.masked_where( xmsk==0, xlon )
    xlat = np.ma.masked_where( xmsk==0, xlat )
    xtim = np.ma.masked_where( xmsk==0, xtim )
        
    (idx,) = np.where(IDs2Follow==id_longest_rec)
    idlngst = idx[0]
        
    for jr in range(nRmax):
        print('')

        # Here records do not correspond to the same time across the buoys!
        # First record of buoys N can be in april, why that of buoy M in january!!!
        
        crec = 'rec%3.3i'%(jr)
        
        rtime = rdt1
        ctime = ''

        Nalive = np.sum(xmsk[:,jr])

        print('\n\n *** Number of buoys still alive at record #'+str(jr)+' => ',Nalive)

        # Trick we put the date vector at `pvIDs=` so the date is writen next to the buoy instead of ID!
        vtim = xtim[:,jr].data ; # `vtim` will contain `0s` rather than '--' for masked records        
        cdates = np.array( [ '('+str(IDs2Follow[i])+')'+mjt.epoch2clock(vtim[i]) for i in range(nBf) ] , dtype='U32' )
        
        ik = mjt.ShowBuoysMap( rtime, xlon[:,jr], xlat[:,jr], pvIDs=cdates, cfig='buoys_RGPS_'+crec+'_'+ctime+'.png',
                               ms=15, ralpha=0.5, lShowDate=False )
    
