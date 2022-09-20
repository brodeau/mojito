#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import argv, exit
from os import path, mkdir
import numpy as nmp
from re import split

from netCDF4 import Dataset

import climporn as cp

import gonzag as gz



idebug = 0

rd_tol = 2. # tolerance distance in km to conclude it's the same buoy

# What to expect in input netCDF file:
list_expected_dim = [ 'time', 'id_buoy' ]
list_expected_var = [ 'time', 'id_buoy', 'latitude', 'longitude' ]
ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!


# Figure stuff:
l_show_IDs_fig = False  ; # annotate ID beside marker in plot...
color_top = 'w'
clr_yellow = '#ffed00'
pt_sz_track = 5
fig_type='png'
rDPI = 150
rzoom = 1.
vfig_size = [ 7.54*rzoom, 7.2*rzoom ]
vsporg = [0., 0., 1., 1.]
col_bg = '#041a4d'
#
# Projection:
#vp =  ['Arctic', 'stere', -60., 40., 122., 57.,    75.,  -12., 10., 'h' ]  # Big Arctic + Northern Atlantic
vp =  ['Arctic', 'stere', -80., 68., 138.5, 62.,    90.,  -12., 10., 'h' ]  # North Pole Arctic (zoom)

    


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


    vIDs = id_in.variables['id_buoy'][:]


    vid_mask = nmp.zeros(Nb, dtype=int)
    vid_mask[:] = 1

    
    # First time step:
    jt   = 0
    vlat = id_in.variables['latitude' ][jt,:]
    vlon = id_in.variables['longitude'][jt,:]

    
    for jb in range(Nb):
        # All the following work should be done if present buoy has not been cancelled yet
        if vid_mask[jb] == 1:
        
            jid = vIDs[jb]
            if idebug>0: print('\n *** DOUBLON scanning: Buoy #'+str(jb)+'=> ID ='+str(jid))
    
            rlat = vlat[jb]
            rlon = vlon[jb]
    
            if idebug>0: print('    ==> lat,lon = ',rlat,rlon)
            
            # Building array of distances (km) with all other buoy at this same time record:
            vdist = gz.Haversine( rlat, rlon, vlat, vlon )
    
            # Need to mask itself (distance = 0!)
            if vdist[jb] != 0.:
                print(' PROBLEM: distance with yourself should be 0!')
            vdist[jb] = 9999.
    
            if idebug>1:
                for rd in vdist:
                    if rd<10.: print('   distance = ', rd)
    
            rdmin = nmp.min(vdist)
            if idebug>0: print('    ==> closest other buoy is at '+str(round(rdmin,2))+' km !')
            
            if rdmin < rd_tol:
                # There is at least 1 buoy too close!
                #  => populating the buoys closer than `rd_tol` km:
                (vidx_close,) = nmp.where( vdist < rd_tol )
    
                if idebug>0:
                    print('    ==> list of buoys closer than '+str(rd_tol)+' km:')
                    for ii in vidx_close:
                        print(vIDs[ii])
    
                vid_mask[(vidx_close,)] = 0



    print('\n LOLO final!')

    # Alright, all IDs spotted via vid_mask should be cancelled in the 2D fields, at all time steps
    xlat = id_in.variables['latitude' ][:,:]
    xlon = id_in.variables['longitude'][:,:]


    vIDs = nmp.ma.masked_where( vid_mask==0, vIDs )
    for jt in range(Nt):
        xlat[jt,:] =  nmp.ma.masked_where( vid_mask==0, xlat[jt,:] )
        xlon[jt,:] =  nmp.ma.masked_where( vid_mask==0, xlon[jt,:] )

            
        
    id_in.close()

