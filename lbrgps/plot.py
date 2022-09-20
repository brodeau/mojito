#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import exit
from os import path, mkdir
import numpy as nmp
from re import split

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import shiftgrid

import climporn as cp

#idebug = 0


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



# Projection:
#vp =  ['Arctic', 'stere', -60., 40., 122., 57.,    75.,  -12., 10., 'h' ]  # Big Arctic + Northern Atlantic
vp =  ['Arctic', 'stere', -80., 68., 138.5, 62.,    90.,  -12., 10., 'h' ]  # North Pole Arctic (zoom)


def ShowBuoysMap( pvt, pvlon, pvlat, pvIDs=[], cnmfig='buoys_RGPS' ):
    '''
        IN:
            * pvt   => vector of length Nt containing the dates as epoch/unix time (integer)
            * pvlon => 2D array (Nt,Nb) of longitudes (float)
            * pvlat => 2D array (Nt,Nb) of  latitudes (float)
            * pvIDs => (OPTIONAL) vector of length Nb of buoys IDs (integer)
    '''
    
    (Nt,Nb) = nmp.shape(pvlon)
    if Nt != len(pvt):
        print("\n *** ERROR [ShowBuoysMap]: record length different for `pvt` and `coordinates`!")
        exit(0)
    
    if not path.exists("./figs"): mkdir("./figs")
    
    cp.fig_style( rzoom, clr_top=color_top )

    PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                   resolution=vp[9], area_thresh=1000., projection='stere', \
                   lat_0=vp[6], lon_0=vp[7], epsg=None)

    Nt = len(pvt)

    for jt in range(Nt):

        ct = cp.epoch2clock(pvt[jt])
    
        print('\n *** [ShowBuoysMap] plotting for time = '+ct)    
        #cfig = 'buoys_'+'%3.3i'%(jt+1)+'.'+fig_type #
        cfig = './figs/'+cnmfig+'_'+split('_',ct)[0]+'.png'

        fig = plt.figure(num=1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
        ax  = plt.axes(vsporg, facecolor=col_bg)

        x0,y0 = PROJ(pvlon[jt,:],pvlat[jt,:])
        #csct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=pt_sz_track*rzoom )
        csct = plt.scatter(x0, y0, marker='o', facecolors='w', edgecolors='none', alpha=0.5, s=pt_sz_track*rzoom ) ; # facecolors='none', edgecolors='r'

        if len(pvIDs) > 0:
            for ii in range(Nb):
                #x0,y0 = PROJ(pvlon[jt,ii],pvlat[jt,ii]) # 
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

    #
    return 0
        
    

