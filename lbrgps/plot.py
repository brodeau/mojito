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
#pt_sz_track = 5
fig_type='png'
rDPI = 150
rzoom = 1.
vfig_size = [ 7.54*rzoom, 7.2*rzoom ]
vsporg = [0., 0., 1., 1.]
col_bg = '#041a4d'

col_red = '#873520'
col_blu = '#3475a3'


msPoints = 20  ; # size  of markers for points aka vertices...
clPoints = 'w' ; # color   "                  "
clPNames = '0.4' ; # color for city/point annotations

# Projection:
#vp =  ['Arctic', 'stere', -60., 40., 122., 57.,    75.,  -12., 10., 'h' ]  # Big Arctic + Northern Atlantic
vp =  ['Arctic', 'stere', -80., 68., 138.5, 62.,    90.,  -12., 10., 'h' ]  # North Pole Arctic (zoom)



def __initStyle__( font_rat=1., color_top='k' ):
    #
    global cfont_clb, cfont_clock, cfont_axis, cfont_ttl, cfont_mail
    #
    params = { 'font.family':'Open Sans',
               'font.weight':    'normal',
               'font.size':       int(15.*font_rat),
               'legend.fontsize': int(22.*font_rat),
               'xtick.labelsize': int(15.*font_rat),
               'ytick.labelsize': int(15.*font_rat),
               'axes.labelsize':  int(17.*font_rat) }
    #
    mpl.rcParams.update(params)
    #
    cfont_clb   = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
    cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(18.*font_rat), 'color':color_top }
    cfont_axis  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
    cfont_ttl   = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(25.*font_rat), 'color':color_top }
    cfont_mail  = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
    #
    return 0





def __figMap__( pt, pvlon, pvlat, BMProj, cdate='', pvIDs=[], cfig='buoys_RGPS.png', ms=5, ralpha=0.5, caller='unknown' ):
    '''
        IN:
            * pt     => the date as epoch/unix time (integer)
            * pvlon  => 1D array (Nb) of longitudes (float)
            * pvlat  => 1D array (Nb) of  latitudes (float)
            * BMProj => the BaseMap projection object
            * cdate  => string of current human-understandable date!
            * pvIDs  => (OPTIONAL) vector of length Nb of buoys IDs (integer)
    '''
    #
    fig = plt.figure(num=1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
    ax  = plt.axes(vsporg, facecolor=col_bg)

    x0,y0 = BMProj(pvlon[:],pvlat[:])
    #csct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=ms*rzoom )
    csct = plt.scatter(x0, y0, marker='o', facecolors='w', edgecolors='none', alpha=ralpha, s=ms*rzoom ) ; # facecolors='none', edgecolors='r'

    # Add IDs figure right next to buoys:
    if len(pvIDs) > 0:
        if Nb != len(pvIDs):
            print('\n *** ERROR ['+caller+'/__figMap__]: `Nb` different for `pvIDs` and `coordinates`!'); exit(0)
        for ii in range(Nb):
            x0,y0 = BMProj(pvlon[jt,ii],pvlat[jt,ii])
            ax.annotate(str(pvIDs[ii]), xy=(x0,y0), xycoords='data', **cp.fig_style.cfont_mrkr)



    BMProj.drawcoastlines(linewidth=0.5)
    BMProj.fillcontinents(color='grey') #, alpha=0)
    #BMProj.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    #BMProj.drawmapboundary()
    BMProj.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1], linewidth=0.3)
    BMProj.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0], linewidth=0.3)

    if cdate != '':
        ax.annotate('Date: '+cdate, xy=(0.6, 0.93), xycoords='figure fraction', **cp.fig_style.cfont_clck)

    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)

    return 0


def ShowBuoysMap( pt, pvlon, pvlat, pvIDs=[], cnmfig='buoys_RGPS', ms=5, ralpha=0.5 ):
    '''
        IN:
            * pt    => the date as epoch/unix time (integer)
            * pvlon => 1D array (Nb) of longitudes (float)
            * pvlat => 1D array (Nb) of  latitudes (float)
            * pvIDs => (OPTIONAL) vector of length Nb of buoys IDs (integer)
    '''
    Nb = len(pvlon)

    if len(pvlat) != Nb:
        print('\n *** ERROR [ShowBuoysMap]: lon and lat vectors disagree in length!')
        exit(0)

    if not path.exists('./figs'): mkdir('./figs')

    cp.fig_style( rzoom, clr_top=color_top )

    PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                   resolution=vp[9], area_thresh=1000., projection='stere', \
                   lat_0=vp[6], lon_0=vp[7], epsg=None)

    ct = cp.epoch2clock(pt)

    print('\n *** [ShowBuoysMap] plotting for time = '+ct)
    cfig = './figs/'+cnmfig+'_'+split('_',ct)[0]+'.png' ; #cfig = 'buoys_'+'%3.3i'%(jt+1)+'.'+fig_type #

    return __figMap__( pt, pvlon, pvlat, PROJ, cdate=ct, pvIDs=pvIDs, cfig=cfig, ms=ms, ralpha=ralpha, caller='ShowBuoysMap' )


def ShowBuoysMap_Trec( pvt, pvlon, pvlat, pvIDs=[], cnmfig='buoys_RGPS', ms=5, ralpha=0.5, clock_res='s' ):
    '''
        IN:
            * pvt   => vector of length Nt containing the dates as epoch/unix time (integer)
            * pvlon => 2D array (Nt,Nb) of longitudes (float)
            * pvlat => 2D array (Nt,Nb) of  latitudes (float)
            * pvIDs => (OPTIONAL) vector of length Nb of buoys IDs (integer)
    '''

    (Nt,Nb) = nmp.shape(pvlon)
    if Nt != len(pvt):
        print('\n *** ERROR [ShowBuoysMap_Trec]: record length different for `pvt` and `coordinates`!')
        exit(0)

    if not path.exists('./figs'): mkdir('./figs')

    cp.fig_style( rzoom, clr_top=color_top )

    PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                   resolution=vp[9], area_thresh=1000., projection='stere', \
                   lat_0=vp[6], lon_0=vp[7], epsg=None)

    Nt = len(pvt)

    kk = 0

    for jt in range(Nt):

        ct = cp.epoch2clock(pvt[jt])

        cfig = './figs/'+cnmfig+'_'+split('_',ct)[0]+'.png' ; #cfig = 'buoys_'+'%3.3i'%(jt+1)+'.'+fig_type #
        
        if clock_res=='d':
            # Daily precision for clock on figure:
            ct = split('_',ct)[0]

        print('\n *** [ShowBuoysMap_Trec] plotting for time = '+ct)

        

        kk = kk + __figMap__( pvt[jt], pvlon[jt,:], pvlat[jt,:], PROJ, cdate=ct, pvIDs=pvIDs, cfig=cfig, ms=ms, ralpha=ralpha )
    #
    return kk



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




def ShowTQMesh( pX, pY, cfig='mesh_quad_map.png', pnames=[], TriMesh=[], QuadMesh=[], lProj=True, zoom=1 ):
    '''
    ### Show points, triangle, and quad meshes on the map!
    ###
    ###  * lProj: True   => we expect degrees for `pX,pY` and use a projection
    ###           False  => we expect km or m for `pX,pY` => cartesian coordinates !
    '''
    from math import log
    
    zrat = 1./zoom
    zrat_log = 1./log(zoom)
    
    kk = __initStyle__(font_rat=zrat_log)    
    
    if lProj:
        # Geograhic coordinates (lon,lat)
        import cartopy.crs as ccrs
        Proj = ccrs.PlateCarree()
        vfig = (12*zoom,9*zoom)
    else:
        # Cartesian coordinates (x,y)
        #  => we want to preserve aspect ratio!
        xA, yA = nmp.min(pX), nmp.min(pY)
        xB, yB = nmp.max(pX), nmp.max(pY)
        Lx, Ly = xB-xA, yB-yA
        dx, dy = 0.05*(Lx), 0.05*(Ly)
        vfig = (10*zoom,10*Ly/Lx*zoom)
    
    fig = plt.figure(num=1, figsize=vfig, facecolor='white')
    
    if lProj:         
        ax   = plt.axes([0.02, 0.02, 0.96, 0.96], projection=Proj)
        ax.stock_img()
        ax.set_extent([-15, 30, 32, 65], crs=Proj) ; #fixme
    else:
        ddx = dx*Ly/Lx
        ax   = plt.axes([ddx/Lx, dy/Ly, (Lx-2*ddx)/Lx, (Ly-2*dy)/Ly], facecolor='0.75')        
        plt.axis([ xA-dx,xB+dx , yA-dy,yB+dy ])

    # Showing points:
    plt.plot( pX, pY, '.', ms=msPoints*zrat, color=clPoints, zorder=200) ; #, alpha=0.5)

    # Adding triangles:
    if len(TriMesh)>0:
        (nbT,_) = nmp.shape(TriMesh); # Number of triangles
        plt.triplot(pX, pY, TriMesh, color=col_red, linestyle='-', lw=2*zrat, zorder=50)    
        # Indicate triangle # in its center:
        for jT in range(nbT):
            vids  = TriMesh[jT,:] ; # the IDs of the 3 points that constitute our triangle        
            rmLon, rmLat = nmp.mean( pX[vids] ), nmp.mean( pY[vids] ) ; # Lon,Lat at center of triangle
            ax.annotate(str(jT), (rmLon, rmLat), color=col_red, fontweight='normal', zorder=60)
    
    # Adding quadrangles:
    if len(QuadMesh)>0:
        (nbQ,_) = nmp.shape(QuadMesh)
        for jQ in range(nbQ):
            vids = QuadMesh[jQ,:]
            vx, vy  = pX[vids], pY[vids]
            vX, vY = nmp.concatenate([vx[:], vx[0:1]]), nmp.concatenate([vy[:],vy[0:1]])  ; # need to make an overlap to close the line
            plt.plot(vX,vY, col_blu, lw=5*zrat, zorder=100)
            plt.fill_between(vX, vY, fc=col_blu, zorder=150, alpha=0.4)
            # Indicate quadrangle # in its center:
            rmLon, rmLat = nmp.mean( pX[vids] ), nmp.mean( pY[vids] ) ; # Lon,Lat at center of triangle
            ax.annotate(str(jQ), (rmLon, rmLat), color='w', fontweight='bold', zorder=160)

    if len(pnames)>0:
        i=0
        for city in pnames:
            ax.annotate(city, (pX[i], pY[i]), color=clPNames, fontweight='bold', zorder=500)
            i=i+1

    plt.savefig(cfig)
    plt.close(1)

