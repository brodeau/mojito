#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import argv, exit
#from os import path
import numpy as nmp
#from re import split

from netCDF4 import Dataset

import climporn as cp
import lbrgps   as lbr

idebug = 0

# What to expect in input netCDF file:
list_expected_dim = [ 'time', 'id_buoy' ]
list_expected_var = [ 'time', 'id_buoy', 'latitude', 'longitude' ]
ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!



if __name__ == '__main__':

    narg = len(argv)
    if not narg in [2]:
        print('Usage: '+argv[0]+' <file_RGPS_common_time.nc>')
        exit(0)
    cf_in  = argv[1]


    # Opening and inspecting input file
    cp.chck4f(cf_in)

    id_in    = Dataset(cf_in)
    #
    list_dim = list( id_in.dimensions.keys() ) ; #print(' ==> dimensions: ', list_dim, '\n')
    for cd in list_expected_dim:
        if not cd in list_dim:
            print(' ERROR: no dimensions `'+cd+'` found into input file!'); exit(0)
    #
    list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
    for cv in list_expected_var:
        if not cv in list_var:
            print(' ERROR: no variable `'+cv+'` found into input file!'); exit(0)
    #
    Nt = id_in.dimensions['time'].size
    print('\n *** Number of records: '+str(Nt))
    Nb = id_in.dimensions['id_buoy'].size
    print('\n *** Number of buoys: '+str(Nb))


    ctunits = id_in.variables['time'].units
    if not ctunits == ctunits_expected:
        print(" ERROR: we expect '"+ctunits_expected+"' as units for time variables, yet we have: "+ctunits)

    vtime = nmp.zeros(Nt, dtype=int)
    vtime[:] = id_in.variables['time'][:]

    cdt1 = cp.epoch2clock(vtime[0])
    cdt2 = cp.epoch2clock(vtime[Nt-1])
    print("\n *** Time range:")
    print(" ==> "+cdt1+" to "+cdt2 )

    # Position of each buoy at this particular time record:
    xlon = id_in.variables['longitude'][:,:]
    xlat = id_in.variables['latitude' ][:,:]

    id_in.close()

    
    kf = lbr.ShowBuoysMap_Trec( vtime, xlon, xlat, ms=5, ralpha=0.5 )
    
