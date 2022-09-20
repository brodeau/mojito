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

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import shiftgrid

import climporn as cp




idebug = 0

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



    # Time for figure:
    ##################

    if not path.exists("./figs"): mkdir("./figs")
        
    cp.fig_style( rzoom, clr_top=color_top )
    
    PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                   resolution=vp[9], area_thresh=1000., projection='stere', \
                   lat_0=vp[6], lon_0=vp[7], epsg=None)

    for jt in range(Nt):

        ct = cp.epoch2clock(vtime[jt])
    
        print('\n *** Plotting for time = '+ct)

        # Position of each buoy at this particular time record:
        vlon = id_in.variables['longitude'][jt,:]
        vlat = id_in.variables['latitude' ][jt,:]
        
        #cfig = 'buoys_'+'%3.3i'%(jt+1)+'.'+fig_type #
        cfig = './figs/buoys_RGPS_'+split('_',ct)[0]+'.png'

        fig = plt.figure(num=1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
        ax  = plt.axes(vsporg, facecolor = col_bg)

        x0,y0 = PROJ(vlon,vlat)
        #csct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=pt_sz_track*rzoom )
        csct = plt.scatter(x0, y0, marker='o', facecolors='w', edgecolors='none', alpha=0.5, s=pt_sz_track*rzoom ) ; # facecolors='none', edgecolors='r'

        if l_show_IDs_fig:
            for ii in range(Nb):
                #if vmask[ii] == 1:
                x0,y0 = PROJ(xlon[jt,ii],xlat[jt,ii]) # 
                ax.annotate(str(vIDs[ii]), xy=(x0,y0), xycoords='data', **cp.fig_style.cfont_mrkr)
        
        PROJ.drawcoastlines(linewidth=0.5)
        PROJ.fillcontinents(color='grey') #, alpha=0)
        #PROJ.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
        #PROJ.drawmapboundary()
        PROJ.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1], linewidth=0.3)
        PROJ.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0], linewidth=0.3)

        ax.annotate('Date: '+ct, xy=(0.6, 0.93), xycoords='figure fraction', **cp.fig_style.cfont_clck)
    
        print('     ===> saving figure: '+cfig)
        plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
        plt.close(1)

        
    id_in.close()

