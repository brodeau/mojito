#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import exit
from os import path, mkdir
import numpy as np
from re import split

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#
from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import shiftgrid

import climporn as cp

#idebug = 0

rndaxiskm = 500 ;

# Figure stuff:
l_show_IDs_fig = False  ; # annotate ID beside marker in plot...
color_top = 'w'
clr_yellow = '#ffed00'

fig_type='png'
rDPI = 150
rzoom = 1.
vfig_size = [ 7.54*rzoom, 7.2*rzoom ]
vsporg = [0., 0., 1., 1.]
col_bg = '#041a4d'

col_red = '#873520'
col_blu = '#3475a3'


msPoints = 8  ; # size  of markers for points aka vertices...
clPoints = 'w' ; # color   "                  "
clPNames = 'w' ; # color for city/point annotations

# Projection:
#vp =  ['Arctic', 'stere', -60., 40., 122., 57.,    75.,  -12., 10., 'h' ]  # Big Arctic + Northern Atlantic
vp =  ['Arctic', 'stere', -80., 68., 138.5, 62.,    90.,  -12., 10., 'h' ]  # North Pole Arctic (zoom)



def _initStyle_( fntzoom=1., color_top='k' ):
    #
    global cfont_clb, cfont_clock, cfont_axis, cfont_ttl, cfont_mail
    #
    fntzoom_inv = 1.*fntzoom**0.5
    #
    params = { 'font.family':'Open Sans',
               'font.weight':    'normal',
               'font.size':       int(18.*fntzoom_inv),
               'legend.fontsize': int(16.*fntzoom_inv),
               'xtick.labelsize': int(15.*fntzoom_inv),
               'ytick.labelsize': int(15.*fntzoom_inv),
               'axes.labelsize':  int(17.*fntzoom) }
    #
    mpl.rcParams.update(params)
    #
    cfont_clb   = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(12.*fntzoom), 'color':color_top }
    cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(16.*fntzoom_inv), 'color':color_top }
    cfont_axis  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*fntzoom), 'color':color_top }
    cfont_ttl   = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(20.*fntzoom), 'color':color_top }
    cfont_mail  = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*fntzoom), 'color':'0.8'}
    #
    return 0

def roundAxisRange( pR, rndKM=None ):
    from math import floor,ceil
    if rndKM:
        zrnd = rndKM
    else:
        zrnd = rndaxiskm
    return [ floor(np.min(pR)/zrnd)*zrnd , ceil(np.max(pR)/zrnd)*zrnd ]

def _set_fig_axis_( pX, pY, rdr=0.05, zoom=1., rangeX=None, rangeY=None ):
    #  => we want to preserve aspect ratio!
    # Rounding axis boundaries at 100 km:
    lFrxy = False
    if rangeX and rangeY: lFrxy = ( len(rangeX)==2 and len(rangeY)==2 )
    if lFrxy:
        xA, xB = rangeX[0], rangeX[1]
        yA, yB = rangeY[0], rangeY[1]
    else:
        [xA, xB] = roundAxisRange(pX)
        [yA, yB] = roundAxisRange(pY)
    #
    Lx, Ly = xB-xA, yB-yA
    #dx, dy = rdr*(Lx)*1.5/zoom, rdr*(Ly)*1.5/zoom
    dx, dy = rdr*Lx, rdr*Ly
    vfig   = (10*zoom,10*Ly/Lx*zoom)
    if not lFrxy:
        xA, xB = xA-dx, xB+dx
        yA, yB = yA-dy, yB+dy
    #
    return (xA,xB), (yA,yB), (Lx,Ly), (dx,dy), vfig




def _figMap_( pt, pvlon, pvlat, BMProj, cdate='', pvIDs=[], cfig='buoys_RGPS.png',
              ms=5, ralpha=0.5, caller='unknown', rzoom=1., title=None ):
    '''
        IN:
            * pt     => the date as epoch/unix time (integer)
            * pvlon  => 1D array (Nb) of longitudes (float)
            * pvlat  => 1D array (Nb) of  latitudes (float)
            * BMProj => the BaseMap projection object
            * cdate  => string of current human-understandable date!
            * pvIDs  => (OPTIONAL) vector of length Nb of buoys IDs (integer)
    '''
    (Nb,) = np.shape(pvlon)
    #
    lIDs = ( np.shape(pvIDs)==(Nb,) )
    #
    # Number of remaining valid points:
    #NbValid = None
    if np.ma.isMaskedArray(pvlat):
        NbValid = pvlon.count()
    else:
        NbValid = 0
        #NbValid = ( (pvlon>=-180.) and (pvlon<=360.)(pvlat>=-90.) and (pvlat<=-90.) ).sum()
    #
    fig = plt.figure(num=1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
    ax  = plt.axes(vsporg, facecolor=col_bg)

    x0,y0 = BMProj(pvlon[:],pvlat[:])
    #csct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=ms*rzoom )
    csct = plt.scatter(x0, y0, marker='o', facecolors='w', edgecolors='none', alpha=ralpha, s=ms*rzoom ) ; # facecolors='none', edgecolors='r'

    # Add IDs figure right next to buoys (if pvIDs provided and not too many buoys!):
    if lIDs and Nb<=200:
        ctype = str(pvIDs.dtype) ; # type of pvIDs:
        lstr = ( ctype[0:2] == '<U' ) ; # does pvIDs contain strings ???
        #
        for ii in range(Nb):
            x0,y0 = BMProj(pvlon[ii],pvlat[ii])
            if lstr:
                ax.annotate(    pvIDs[ii] , xy=(x0,y0), xycoords='data', **cp.fig_style.cfont_mrkr)
            else:
                ax.annotate(str(pvIDs[ii]), xy=(x0,y0), xycoords='data', **cp.fig_style.cfont_mrkr)



    BMProj.drawcoastlines(linewidth=0.5)
    BMProj.fillcontinents(color='grey') #, alpha=0)
    #BMProj.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    #BMProj.drawmapboundary()
    BMProj.drawmeridians(np.arange(-180,180,20), labels=[0,0,0,1], linewidth=0.3)
    BMProj.drawparallels(np.arange( -90, 90,10), labels=[1,0,0,0], linewidth=0.3)

    if cdate != '':
        ax.annotate('Date: '+cdate, xy=(0.6, 0.925), xycoords='figure fraction', **cp.fig_style.cfont_clck)

    if NbValid:
        ax.annotate(' Nb. buoys = '+str(NbValid), xy=(0.02, 0.8), xycoords='figure fraction', **cp.fig_style.cfont_clck)

    if title:
        ax.annotate(title, xy=(0.5, 0.965), xycoords='figure fraction', ha='center', **cp.fig_style.cfont_ttl)

    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)

    return 0


def ShowBuoysMap( pt, pvlon, pvlat, pvIDs=[], cfig='buoys_RGPS.png', cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1., title=None ):
    '''
        IN:
            * pt    => the date as epoch/unix time (integer)
            * pvlon => 1D array (Nb) of longitudes (float)
            * pvlat => 1D array (Nb) of  latitudes (float)
            * pvIDs => (OPTIONAL) vector of length Nb of buoys IDs (integer)
    '''
    (Nb,) = np.shape(pvlon)
    if np.shape(pvlat) != (Nb,):
        print('\n *** ERROR [ShowBuoysMap]: lon and lat vectors disagree in shape!'); exit(0)

    if not path.exists('./figs'): mkdir('./figs')

    cp.fig_style( zoom, clr_top=color_top )

    PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                   resolution=vp[9], area_thresh=1000., projection='stere', \
                   lat_0=vp[6], lon_0=vp[7], epsg=None)

    if lShowDate:
        ct = cp.epoch2clock(pt)
    else:
        ct = ''

    if lShowDate:  print('\n *** [ShowBuoysMap] plotting for time = '+ct)
    if cnmfig:
        cfig = './figs/'+cnmfig+'_'+split('_',ct)[0]+'.png' ; #cfig = 'buoys_'+'%3.3i'%(jt+1)+'.'+fig_type #

    return _figMap_( pt, pvlon, pvlat, PROJ, cdate=ct, pvIDs=pvIDs, cfig=cfig, ms=ms,
                     ralpha=ralpha, caller='ShowBuoysMap', rzoom=zoom, title=title )



def ShowBuoysMap_Trec( pvt, pvlon, pvlat, pvIDs=[], cnmfig='buoys_RGPS', ms=5, ralpha=0.5, clock_res='s', NminPnts=100 ):
    '''
        IN:
            * pvt   => vector of length Nt containing the dates as epoch/unix time (integer)
            * pvlon => 2D array (Nt,Nb) of longitudes (float)
            * pvlat => 2D array (Nt,Nb) of  latitudes (float)
            * pvIDs => (OPTIONAL) vector of length Nb of buoys IDs (integer)
            *
            * NminPnts => if a record contains less points than this number we stop plotting!
    '''

    (Nt,Nb) = np.shape(pvlon)
    if Nt != len(pvt):
        print('\n *** ERROR [ShowBuoysMap_Trec]: record length different for `pvt` and `coordinates`!',len(pvt),Nt)
        exit(0)

    if not path.exists('./figs'): mkdir('./figs')

    cp.fig_style( rzoom, clr_top=color_top )

    PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                   resolution=vp[9], area_thresh=1000., projection='stere', \
                   lat_0=vp[6], lon_0=vp[7], epsg=None)

    Nt = len(pvt)

    kk = 0

    for jt in range(Nt):

        #  How many buoys alive at this record ?
        (idx,) = np.where( pvlat[jt,:] >= -90. ) ; # ( since `pvlat` should be masked with `-9999` values...)
        nPtsAlive = len(idx)

        if nPtsAlive >= NminPnts:

            ct = cp.epoch2clock(pvt[jt])

            cfig = './figs/'+cnmfig+'_'+split('_',ct)[0]+'.png' ; #cfig = 'buoys_'+'%3.3i'%(jt+1)+'.'+fig_type #

            if clock_res=='d':
                # Daily precision for clock on figure:
                ct = split('_',ct)[0]

            print('\n *** [ShowBuoysMap_Trec] plotting for time = '+ct)

            kk = kk + _figMap_( pvt[jt], pvlon[jt,:], pvlat[jt,:], PROJ, cdate=ct, pvIDs=pvIDs, cfig=cfig, ms=ms, ralpha=ralpha )
    #
    return kk



def plot_interp_series( iID, cname, vTs, vTt, vFs, vFt ):
    #
    # For debugging:
    # Overlay of original source F(t) series vFs(vTs)
    # and interpolated version vFt on vTt axis
    #
    import matplotlib.dates as mdates
    from datetime import datetime as dtm
    #
    cfig = 'debug_interp_'+cname+'_ID'+'%6.6i'%(iID)+'.'+fig_type
    fig = plt.figure(num=1, figsize=(12,5), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes([0.05, 0.13, 0.9, 0.8])
    #lolo
    #
    #
    func = np.vectorize(dtm.utcfromtimestamp)
    ax.set_xlim( dtm.utcfromtimestamp(vTt[0,1]), dtm.utcfromtimestamp(vTt[-1,2]) )
    pl1 = plt.plot( mdates.date2num(func(vTs)), vFs, 'o-', color='#041a4d', linewidth=4., ms=10.)
    pl2 = plt.plot( mdates.date2num(func(vTt[:,0])), vFt, '*-', color='r'      , linewidth=1., ms=3)
    date_fmt = '%Y/%m/%d'
    date_formatter = mdates.DateFormatter(date_fmt)
    ax.xaxis.set_major_formatter(date_formatter)
    fig.autofmt_xdate()
    #
    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)




def ShowTQMesh( pX, pY, cfig='mesh_quad_map.png', pnames=[], ppntIDs=[], qIDs=[], qnames=[], TriMesh=[],
                pX_Q=[], pY_Q=[], QuadMesh=[], lGeoCoor=True, zoom=1, cProj='NPS',
                rangeX=None, rangeY=None ):
    '''
    ### Show points, triangle, and quad meshes on the map!
    ###
    ###  * lGeoCoor: True   => we expect degrees for `pX,pY` => geographic (lon,lat) coordinates !
    ###           False  => we expect km or m for `pX,pY` => cartesian coordinates !
    ###
    ###     Specify pX_Q & pY_Q when plotting QuadMesh when IDs are not those of the
    ###     traingle world!
    ###
    ###  * pnames:  (len=nP) name (string) for each point
    ###  * ppntIDs: (len=nP) ID (integer)  for each point
    '''
    zrat = 1./zoom
    kk = _initStyle_(fntzoom=zoom)

    rsz_annot  = 5*zoom**0.4

    (nbP,) = np.shape(pX) ; # Number of points that defines all the triangles...

    if lGeoCoor:
        # Geograhic coordinates (lon,lat)
        import cartopy.crs     as ccrs
        import cartopy.feature as cftr
        #
        ProjPC = ccrs.PlateCarree()
        #
        vfig = (10*zoom,10*zoom)
        if cProj=='NPS':
            Proj = ccrs.NorthPolarStereo()
        elif cProj=='PC':
            Proj = ProjPC
            vfig = (12*zoom,9*zoom)
        else:
            print('\n *** ERROR ['+caller+'.ShowTQMesh()]: Unknown projection "'+cProj+'"!'); exit(0)
        #
    else:
        # Cartesian coordinates (x,y) in km...
        (xA,xB), (yA,yB), (Lx,Ly), (dx,dy), vfig = _set_fig_axis_( pX, pY, zoom=zoom, rangeX=rangeX, rangeY=rangeY )

    fig = plt.figure(num=1, figsize=vfig, facecolor='w')

    if lGeoCoor:
        ax   = plt.axes([0.02, 0.02, 0.96, 0.96], projection=Proj, facecolor=col_bg)
        ax.set_extent([-180, 180, 65, 90], ProjPC) ; # Alwasy PlateCaree here !!!
        ax.add_feature(cftr.LAND, color='0.5', zorder=70)
        # Showing points:
        plt.plot( pX, pY, '.', ms=msPoints, color=clPoints, zorder=200, transform=ProjPC) ; #, alpha=0.5)  ; # Alwasy PlateCaree here !!!
    else:
        ddx = dx*Ly/Lx
        ax  = plt.axes([1.25*ddx/Lx, 1.25*dy/Ly, (Lx-2*ddx)/Lx, (Ly-2*dy)/Ly], facecolor='0.75')
        plt.axis([ xA,xB , yA,yB ])
        # Showing points:
        plt.plot( pX, pY, '.', ms=msPoints, color=clPoints, zorder=200 )

    # Adding triangles:
    if len(TriMesh)>0:
        (nbT,_) = np.shape(TriMesh); # Number of triangles
        if lGeoCoor:
            plt.triplot(pX, pY, TriMesh, color=col_red, linestyle='-', lw=2*zrat, zorder=50, transform=ProjPC)
        else:
            plt.triplot(pX, pY, TriMesh, color=col_red, linestyle='-', lw=2*zrat, zorder=50)

        # Indicate triangle # in its center:
        for jT in range(nbT):
            vids  = TriMesh[jT,:] ; # the IDs of the 3 points that constitute our triangle
            rmLon, rmLat = np.mean( pX[vids] ), np.mean( pY[vids] ) ; # Lon,Lat at center of triangle
            ax.annotate(str(jT), (rmLon, rmLat), color=col_red, fontweight='normal', size=rsz_annot, zorder=60)

    # Adding quadrangles:
    if len(QuadMesh)>0:
        (nbQ,_) = np.shape(QuadMesh)

        for jQ in range(nbQ):
            vids = QuadMesh[jQ,:] ; # the 4 IDs of the 4 points defining this Quad
            if len(pX_Q)>0 and len(pY_Q)>0:
                vx, vy = pX_Q[vids], pY_Q[vids]
            else:
                vx, vy = pX[vids], pY[vids]
            vX, vY = np.concatenate([vx[:], vx[0:1]]), np.concatenate([vy[:],vy[0:1]])  ; # need to make an overlap to close the line
            plt.plot(vX,vY, col_blu, lw=3*zrat, zorder=100)
            plt.fill_between(vX, vY, fc=col_blu, zorder=150, alpha=0.4)
            if len(qIDs)>0:
                # Indicate quadrangle # in its center:
                rmLon, rmLat = np.mean( vx ), np.mean( vy ) ; # Lon,Lat at center of triangle
                ax.annotate(qIDs[jQ], (rmLon, rmLat), color='w', fontweight='bold', size=rsz_annot, zorder=160)
            if len(qnames)>0:
                rmLon, rmLat = np.mean( vx ), np.mean( vy ) ; # Lon,Lat at center of triangle
                ax.annotate(qnames[jQ], (rmLon, rmLat), color='b', size=0.6*rsz_annot, zorder=165)

    if len(pnames)>0:
        for jP in range(nbP):
            ax.annotate(pnames[jP], (pX[jP], pY[jP]), color=clPNames, fontweight='bold', size=rsz_annot, zorder=500)
    if len(ppntIDs)>0:
        for jP in range(nbP):
            ax.annotate(str(ppntIDs[jP]), (pX[jP], pY[jP]), color='k', fontweight='bold', size=0.75*rsz_annot, zorder=520)

    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig)
    plt.close(1)






def ShowDeformation( pX, pY, pF, cfig='deformation_map.png', cwhat='div', zoom=1,
                     pFmin=-1., pFmax=1., rangeX=None, rangeY=None, unit=None,
                     marker_size=None, title=None ):
    '''
    ### Show points, triangle, and quad meshes on the map!
    ###
    ###  * lGeoCoor: True   => we expect degrees for `pX,pY` => geographic (lon,lat) coordinates !
    ###           False  => we expect km or m for `pX,pY` => cartesian coordinates !
    ###
    ###     Specify pX_Q & pY_Q when plotting QuadMesh when IDs are not those of the
    ###     traingle world!
    '''
    zrat = 1./zoom
    kk = _initStyle_(fntzoom=zoom)

    # Colormap:
    if   cwhat=='shr':
        cm = plt.cm.get_cmap('viridis')
    elif   cwhat=='tot':
        cm = plt.cm.get_cmap('inferno')
    elif   cwhat=='UMc':
        cm = plt.cm.get_cmap('plasma')
    else:
        cm = plt.cm.get_cmap('RdBu')
    cn = colors.Normalize(vmin=pFmin, vmax=pFmax, clip = False)

    # Cartesian coordinates (x,y)
    (xA,xB), (yA,yB), (Lx,Ly), (dx,dy), vfig = _set_fig_axis_( pX, pY, zoom=zoom, rangeX=rangeX, rangeY=rangeY )

    fig = plt.figure(num=1, figsize=vfig, facecolor='white')

    zix = 1.5*dx/Lx
    ziy = 1.5*dy/Ly
    if unit: ziy = 4.2*dy/Ly
    ax = plt.axes( [ zix, ziy, 1.-(dx/Lx+zix), 1.-(dy/Ly+ziy) ], facecolor='0.75')

    plt.axis([ xA,xB , yA,yB ])

    # Pixel size:
    if not marker_size:
        marker_size = 500.
        if rangeX and rangeY:
            rrm = int( max( abs(rangeX[1]-rangeX[0]) , abs(rangeY[1]-rangeY[0]) ) )
            marker_size = 1710000./rrm
    marker_size = marker_size*500./Ly

    # Showing points:
    plt.scatter( pX, pY, c=pF, s=marker_size, marker='s', cmap=cm, norm=cn )

    if title:
        ax.annotate(title, xy=(0.1, ziy+0.05), xycoords='figure fraction', **cfont_ttl) ; #ha='center'
    if unit:
        # => triggers the colorbar
        ax2 = plt.axes([0.1, ziy/2., 0.8, 0.02])
        clb = mpl.colorbar.ColorbarBase(ax=ax2, cmap=cm, norm=cn, orientation='horizontal', extend='both')
        clb.set_label(unit, **cfont_clb)

    plt.savefig(cfig)
    plt.close(1)
    return 0




def Plot1Mesh( pcoor_trg, Ys, Xs, isrc_msh, vnames=['P1','P2','P3','P4'], pcoor_cntr=(),
               fig_name='mesh.png', pcoor_extra=(-999.,-999.), label_extra=None ):
    '''
    isrc_msh: 2D integer array of shape (4,2)
    pcoor_cntr : (lat,lon) of T (or F) point at center of F-defined (or T-defined) cell
    wghts:    1D real array of shape (4,)
    pcoor_extra: just an extra point to show on the figure....
    '''
    (yT,xT)                             = pcoor_trg
    [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = isrc_msh[:,:]    
    #
    fig = plt.figure(num = 1, figsize=[14,10], facecolor='w', edgecolor='k')
    ax = plt.axes([0.09, 0.07, 0.6, 0.9])

    plt.plot( [     xT    ], [     yT    ], marker='o', ms=15, color='k', label='target point' ) ; # target point !

    for i in range(4):
        [ja,ia] = isrc_msh[i,:]
        [jb,ib] = isrc_msh[(i+1)%4,:]
        # Lines to see the mesh:
        plt.plot( [ Xs[ja,ia],Xs[jb,ib] ], [ Ys[ja,ia],Ys[jb,ib] ], linestyle='-', color='k', marker=None, ms=10, label=None)
        # The point:
        plt.plot( [ Xs[ja,ia] ], [ Ys[ja,ia] ], marker='o', ms=10, label=vnames[i] )
        
    if len(pcoor_cntr)==2:
        (yC,xC) = pcoor_cntr
        plt.plot( [     xC    ], [     yC    ], marker='o', ms=10, color='k', label='center' ) ; # target point !
        
    if pcoor_extra!=(-999.,-999.):
        (yE,xE) = pcoor_extra
        plt.plot( [     xE    ], [     yE    ], marker='+', ms=20, color='0.5', label=label_extra ) ; # target point !

    ax.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), fancybox=True)
    plt.savefig(fig_name, dpi=100, transparent=False)
    plt.close(1)








def PlotPDFdef( pbinb, pbinc, ppdf, Np=None, name='Divergence', cfig='PDF.png',
                xrng=None, title=None, period=None, dx=0.01 ):
    '''
      * pbinb: vector of the bounds of the bins (x-axis), size = nB+1
      * pbinc: vector of the center of the bins (x-axis), size = nB
      * ppdf:  the PDF, same size as `pbinc`, size = nB
    '''
    from math import ceil

    (nB,) = np.shape(ppdf)
    if len(pbinc) != nB:
        print('\n *** ERROR ['+caller+'/PlotPDFdef]: wrong size for `pbinc`!'); exit(0)
    if len(pbinb) != nB+1:
        print('\n *** ERROR ['+caller+'/PlotPDFdef]: wrong size for `pbinb`!'); exit(0)

    xmax = 0.1

    ki = _initStyle_()

    fig = plt.figure( num = 1, figsize=(10,9), dpi=None )
    #
    ax = plt.axes([0.1, 0.085, 0.85, 0.85])
    #
    #ax.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    #
    # Y-axis:
    #rinc = 1. ; # => dy of 0.01
    #ymax = ceil(rinc*np.max(ppdf))/rinc
    #plt.yticks( np.arange(0.,ymax+1./rinc,1./rinc) )
    #
    ymax = 180. ; rinc= 20.
    plt.yticks( np.arange(0.,ymax+rinc,rinc) )
    #
    ax.set_ylim(0.,ymax)
    #ax.set_ylim(0.,np.max(ppdf))
    ax.set_ylabel(r'PDF')
    #
    # X-axis:
    plt.xticks( np.arange(0., xmax+dx, dx) )
    ax.set_xlim( pbinb[0], xmax)
    #ax.set_xlim( 0., xmax)
    plt.xlabel(r''+name+' [day$^{-1}$]', color='k')
    #
    #if nx_subsamp>1:
    #    ax_lab = []
    #    locs, labels = plt.xticks()
    #    cpt = 0 # with ipct = 1: tick priting will start at y1+dt on x axis rather than y1
    #    for rr in locs:
    #        if cpt % nx_subsamp == 0:
    #            if rr%1.0 == 0.:
    #                cr = str(int(rr))  # it's something like 22.0, we want 22 !!!
    #            else:
    #                cr = str(round(rr,6))
    #            ax_lab.append(cr)
    #        else:
    #            ax_lab.append(' ')
    #        cpt = cpt + 1
    #    plt.xticks(locs,ax_lab)
    #    del ax_lab

    #plt.bar ( pbinc[:],  ppdf[:], width=pbinb[1:]-pbinb[:-1], color='0.6', edgecolor=None, linewidth=2, zorder=10 )
    #plt.step( pbinb[1:], ppdf[:],  color='k', linewidth=0.7, zorder=10) ; # the enveloppe
    plt.bar ( pbinc[:],  ppdf[:], width=pbinb[1:]-pbinb[:-1], color='0.6', edgecolor='k', linewidth=0.5, zorder=10 )

    if Np:
        ax.annotate('N = '+str(Np), xy=(0.72, 0.85), xycoords='figure fraction', **cfont_clock)
    #if wbin:
    #    ax.annotate('Bin width = '+str(wbin)+r' day$^{-1}$', xy=(0.62, 0.82), xycoords='figure fraction', **cfont_clock)
    if period:
        ax.annotate('Period = '+period, xy=(0.62, 0.79), xycoords='figure fraction', **cfont_clock)
    if title:
        ax.annotate(title, xy=(0.5, 0.95), xycoords='figure fraction', ha='center', **cfont_ttl)

    plt.savefig(cfig, dpi=100, orientation='portrait', transparent=False)
    plt.close(1)
    print(' * [PlotPDFdef()]: created figure '+cfig)
    return 0




def LogPDFdef( pbinb, pbinc, ppdf, Np=None, name='Divergence', cfig='PDF.png', reskm=None,
               title=None, period=None, origin=None,
               ppdf2=[], Np2=None, origin2=None,
               ppdf3=[], Np3=None, origin3=None ):
    '''
      * pbinb: vector of the bounds of the bins (x-axis), size = nB+1
      * pbinc: vector of the center of the bins (x-axis), size = nB
      * ppdf:  the PDF, same size as `pbinc`, size = nB
    '''
    from math import ceil

    (nB,) = np.shape(ppdf)
    if len(pbinc) != nB:
        print('\n *** ERROR ['+caller+'/LogPDFdef]: wrong size for `pbinc`!'); exit(0)
    if len(pbinb) != nB+1:
        print('\n *** ERROR ['+caller+'/LogPDFdef]: wrong size for `pbinb`!'); exit(0)

    rycut_tiny =1.e-6 ; # mask probability values (Y-axis) below this limit for non-RGPS data!!!

    l_comp2 = ( np.shape(ppdf2)==(nB,) )
    l_comp3 = ( l_comp2 and np.shape(ppdf3)==(nB,) )
    
    # For figure axes:
    xlog_min, xlog_max = 2.75e-3, 0.5
    ylog_min,ylog_max = 5.e-3, 3.5e2
    rxlabs = [0.005, 0.01, 0.05, 0.1, 0.5]
    cxlabs = ['0.005', '0.01', '0.05', '0.1', '0.5']
    lmask_tiny, lmask_tiny2, lmask_tiny3 = (origin!='RGPS'), (origin2!='RGPS'), (origin3!='RGPS')
    
    if reskm:
        if reskm>30 and reskm<70:         
            ylog_min = 3.e-2
            xlog_max = 0.2
            rxlabs = [0.005, 0.01, 0.025, 0.05, 0.1, 0.2]
            cxlabs = ['0.005', '0.01', '0.025', '0.05', '0.1', '0.2']
            rycut_tiny =3.e-2
        if reskm>=70 and reskm<120:         
            ylog_min = 5.e-2
            xlog_max = 0.15
            rxlabs = [0.005, 0.01, 0.025, 0.05, 0.1]
            cxlabs = ['0.005', '0.01', '0.025', '0.05', '0.1']
            rycut_tiny =0.5e-1
        elif reskm>=120 and reskm<240:
            ylog_min = 1.e-1
            xlog_max = 0.1
            rxlabs = [0.005, 0.01, 0.025, 0.05, 0.1]
            cxlabs = ['0.005', '0.01', '0.025', '0.05', '0.1']
            rycut_tiny =1.e-1
        elif reskm>=240 and reskm<480:
            ylog_min = 0.8e-1
            xlog_max = 0.1
            rxlabs = [0.005, 0.01, 0.025, 0.05, 0.1]
            cxlabs = ['0.005', '0.01', '0.025', '0.05', '0.1']
            rycut_tiny =0.7e-1

        if reskm>30:
            lmask_tiny, lmask_tiny2, lmask_tiny3 = True, True, True
    
    ki = _initStyle_()

    fig = plt.figure( num = 1, figsize=(10,9), dpi=None )
    ax = plt.axes([0.11, 0.085, 0.85, 0.85])

    if lmask_tiny: ppdf = np.ma.masked_where( ppdf<rycut_tiny, ppdf )

    clbl = origin
    if Np and origin:
        clbl = origin+' (N = '+str(Np)+')'

    plt.loglog(pbinc[:], ppdf[:], 'o', markersize=12, linestyle='-', linewidth=6, fillstyle='none', color='k', label=clbl, zorder=5)

    if l_comp2:
        clbl = origin2
        if Np2 and origin2:
            origin2 = str.replace( str.replace( origin2, 'NEMO-','') , '_NANUK4', '')
            clbl = origin2+' (N = '+str(Np2)+')'
        if lmask_tiny2: ppdf2 = np.ma.masked_where( ppdf2<rycut_tiny, ppdf2 )
        plt.loglog(pbinc[:], ppdf2[:], 's', markersize=12, fillstyle='none', color='0.4', linestyle='-', linewidth=4,  label=clbl, zorder=10)
        ax.legend(loc='center left', fancybox=True) ; # , bbox_to_anchor=(1.07, 0.5)

    if l_comp3:
        clbl = origin3
        if Np3 and origin3:
            origin3 = str.replace( str.replace( origin3, 'NEMO-','') , '_NANUK4', '')
            clbl = origin3+' (N = '+str(Np3)+')'
        if lmask_tiny3: ppdf3 = np.ma.masked_where( ppdf3<rycut_tiny, ppdf3 )
        plt.loglog(pbinc[:], ppdf3[:], '*', markersize=14, color='0.65', linestyle='-', linewidth=4, label=clbl, zorder=10)
        ax.legend(loc='lower left', fancybox=True) ; # , bbox_to_anchor=(1.07, 0.5)


    # X-axis:
    plt.xlabel(r''+name+' [day$^{-1}$]', color='k')
    ax.set_xlim(xlog_min, xlog_max)
    ax.set_xticks(rxlabs)
    ax.set_xticklabels(cxlabs)
    #plt.tick_params(axis='x', which='minor')
    #from matplotlib.ticker import FormatStrFormatter
    #ax.xaxis.set_minor_formatter(FormatStrFormatter("%.3f"))

    # Y-axis:
    plt.ylabel('PDF', color='k')
    ax.set_ylim(ylog_min, ylog_max)

    ax.grid(color='0.5', linestyle='-', which='minor', linewidth=0.2, zorder=0.1)
    ax.grid(color='0.5', linestyle='-', which='major', linewidth=0.4, zorder=0.1)

    if Np and not (l_comp2 or l_comp3):
        ax.annotate('N = '+str(Np), xy=(0.72, 0.85), xycoords='figure fraction', **cfont_clock)
    if period:
        ax.annotate('Period = '+period, xy=(0.62, 0.79), xycoords='figure fraction', **cfont_clock)
    if title:
        ax.annotate(title, xy=(0.5, 0.95), xycoords='figure fraction', ha='center', **cfont_ttl)

    plt.savefig(cfig, dpi=100, orientation='portrait', transparent=False)
    plt.close(1)
    print(' * [LogPDFdef()]: created figure '+cfig)
    return 0




def ShowDefQuad( pX4, pY4, pF, cfig='deformation_map.png', cwhat='div', zoom=1,
                 pFmin=-1., pFmax=1., rangeX=None, rangeY=None, unit=None, title=None ):
    '''
    ### Show points, triangle, and quad meshes on the map!
    ### => each quadrangle is filled with the appropriate color from colormap !!!
    ###
    ###
    ###  * pX4, pY4: for each quad the coordinates of the 4 vertices!
    ###
    ###  * lGeoCoor: True   => we expect degrees for `pX4,pY4` => geographic (lon,lat) coordinates !
    ###           False  => we expect km or m for `pX4,pY4` => cartesian coordinates !
    ###
    ###     Specify pX4_Q & pY4_Q when plotting QuadMesh when IDs are not those of the
    ###     traingle world!
    '''
    (nQ,) = np.shape(pF)
    if np.shape(pX4)!=(nQ,4) or np.shape(pY4)!=(nQ,4):
        print('\n *** ERROR [ShowDefQuad]: wrong shape for `pX4` or/and `pY4`!'); exit(0)

    kk = _initStyle_(fntzoom=zoom)

    # Colormap:
    if   cwhat=='shr':
        cm = plt.cm.get_cmap('viridis')
    elif   cwhat=='tot':
        cm = plt.cm.get_cmap('inferno')
    elif   cwhat=='UMc':
        cm = plt.cm.get_cmap('plasma')
    else:
        cm = plt.cm.get_cmap('RdBu')
    cn = colors.Normalize(vmin=pFmin, vmax=pFmax, clip = False)

    # Cartesian coordinates (x,y)
    (xA,xB), (yA,yB), (Lx,Ly), (dx,dy), vfig = _set_fig_axis_( np.mean( pX4[:,:], axis=1 ), np.mean( pY4[:,:], axis=1 ),
                                                               zoom=zoom, rangeX=rangeX, rangeY=rangeY )

    fig = plt.figure(num=1, figsize=vfig, facecolor='white')

    zix = 1.5*dx/Lx
    ziy = 1.5*dy/Ly
    if unit: ziy = 4.2*dy/Ly
    ax = plt.axes( [ zix, ziy, 1.-(dx/Lx+zix), 1.-(dy/Ly+ziy) ], facecolor='0.75')

    plt.axis([ xA,xB , yA,yB ])

    for jQ in range(nQ):
        if not np.isnan(pF[jQ]):
            znorm = cn(pF[jQ])
            colrgb = cm(znorm)
            plt.fill( pX4[jQ,:], pY4[jQ,:], facecolor=colrgb, edgecolor=None, linewidth=0. )

    if title:
        ax.annotate(title, xy=(0.1, ziy+0.05), xycoords='figure fraction', **cfont_ttl) ; #ha='center'
    if unit:
        # => triggers the colorbar
        ax2 = plt.axes([0.1, ziy/2., 0.8, 0.02])
        clb = mpl.colorbar.ColorbarBase(ax=ax2, cmap=cm, norm=cn, orientation='horizontal', extend='both')
        clb.set_label(unit, **cfont_clb)

    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig)
    plt.close(1)
    return 0





def ShowDefQuadGeoArctic( pX4, pY4, pF, cfig='deformation_map.png', cwhat='div', zoom=1,
                          pFmin=-1., pFmax=1., rangeX=None, rangeY=None, unit=None, title=None ):
    '''
    ### Show points, triangle, and quad meshes on the map!
    ### => each quadrangle is filled with the appropriate color from colormap !!!
    ###
    ###
    ###  * pX4, pY4: for each quad the coordinates of the 4 vertices!
    ###
    ###  * lGeoCoor: True   => we expect degrees for `pX4,pY4` => geographic (lon,lat) coordinates !
    ###           False  => we expect km or m for `pX4,pY4` => cartesian coordinates !
    ###
    ###     Specify pX4_Q & pY4_Q when plotting QuadMesh when IDs are not those of the
    ###     traingle world!
    '''
    from .util import ConvertCartesianNPSkm2Geo
    #from cartopy.crs import PlateCarree, NorthPolarStereo
    #

    (nQ,) = np.shape(pF)
    if np.shape(pX4)!=(nQ,4) or np.shape(pY4)!=(nQ,4):
        print('\n *** ERROR [ShowDefQuad]: wrong shape for `pX4` or/and `pY4`!'); exit(0)

    kk = _initStyle_(fntzoom=zoom)

    # Need lat,lon from pX4, pY4 !
    #crs_src = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    #crs_trg = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    
    zlat, zlon = ConvertCartesianNPSkm2Geo( pY4, pX4 )

    # Colormap:
    if   cwhat=='shr':
        cm = plt.cm.get_cmap('viridis')
    elif   cwhat=='tot':
        cm = plt.cm.get_cmap('inferno')
    elif   cwhat=='UMc':
        cm = plt.cm.get_cmap('plasma')
    else:
        cm = plt.cm.get_cmap('RdBu')
    cn = colors.Normalize(vmin=pFmin, vmax=pFmax, clip = False)


    PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                   resolution=vp[9], area_thresh=1000., projection='stere', \
                   lat_0=vp[6], lon_0=vp[7], epsg=None)


    fig = plt.figure(num=1, figsize=(7.54*1.3, 7.2*1.3), dpi=None, facecolor=col_bg, edgecolor=col_bg)
    ax  = plt.axes(vsporg, facecolor='w')

    x0,y0 = PROJ(zlon,zlat)

    for jQ in range(nQ):
        if not np.isnan(pF[jQ]):
            znorm = cn(pF[jQ])
            colrgb = cm(znorm)
            plt.fill( x0[jQ,:], y0[jQ,:], facecolor=colrgb, edgecolor=None, linewidth=0. )

    PROJ.drawcoastlines(linewidth=0.5)
    PROJ.fillcontinents(color='grey') #, alpha=0)
    PROJ.drawmeridians(np.arange(-180,180,20), labels=[0,0,0,1], linewidth=0.3)
    PROJ.drawparallels(np.arange( -90, 90,10), labels=[1,0,0,0], linewidth=0.3)


            
    if title:
        ax.annotate(title, xy=(0.6, 0.95), xycoords='figure fraction', **cfont_ttl) ; #ha='center'
    #if unit:
    #    # => triggers the colorbar
    #    ax2 = plt.axes([0.1, ziy/2., 0.8, 0.02])
    #    clb = mpl.colorbar.ColorbarBase(ax=ax2, cmap=cm, norm=cn, orientation='horizontal', extend='both')
    #    clb.set_label(unit, **cfont_clb)

    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig)
    plt.close(1)
    return 0


