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
               'legend.fontsize': int(22.*fntzoom_inv),
               'xtick.labelsize': int(15.*fntzoom_inv),
               'ytick.labelsize': int(15.*fntzoom_inv),
               'axes.labelsize':  int(17.*fntzoom) }
    #
    mpl.rcParams.update(params)
    #
    cfont_clb   = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(12.*fntzoom), 'color':color_top }
    cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(16.*fntzoom_inv), 'color':color_top }
    cfont_axis  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*fntzoom), 'color':color_top }
    cfont_ttl   = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(23.*fntzoom), 'color':color_top }
    cfont_mail  = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*fntzoom), 'color':'0.8'}
    #
    return 0

def _rnd_ax_range_( pR ):
    from math import floor,ceil
    return floor(np.min(pR)/rndaxiskm)*rndaxiskm , ceil(np.max(pR)/rndaxiskm)*rndaxiskm

def _set_fig_axis_( pX, pY, rdr=0.05, zoom=1., rangeX=None, rangeY=None ):    
    #  => we want to preserve aspect ratio!
    # Rounding axis boundaries at 100 km:
    lFrxy = False
    if rangeX and rangeY: lFrxy = ( len(rangeX)==2 and len(rangeY)==2 )
    if lFrxy:
        xA, xB = rangeX[0], rangeX[1]
        yA, yB = rangeY[0], rangeY[1]
    else:
        xA, xB = _rnd_ax_range_(pX)
        yA, yB = _rnd_ax_range_(pY)
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




def _figMap_( pt, pvlon, pvlat, BMProj, cdate='', pvIDs=[], cfig='buoys_RGPS.png', ms=5, ralpha=0.5, caller='unknown', rzoom=1. ):
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
    NbnotMasked = None
    if np.ma.isMaskedArray(pvlon): NbnotMasked = pvlon.count()
    #
    fig = plt.figure(num=1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
    ax  = plt.axes(vsporg, facecolor=col_bg)

    x0,y0 = BMProj(pvlon[:],pvlat[:])
    #csct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=ms*rzoom )
    csct = plt.scatter(x0, y0, marker='o', facecolors='w', edgecolors='none', alpha=ralpha, s=ms*rzoom ) ; # facecolors='none', edgecolors='r'

    # Add IDs figure right next to buoys (if pvIDs provided and not too many buoys!):
    if np.shape(pvIDs)==(Nb,) and Nb<=800:
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
        ax.annotate('Date: '+cdate, xy=(0.6, 0.93), xycoords='figure fraction', **cp.fig_style.cfont_clck)

    if NbnotMasked:
        ax.annotate(' Nb. buoys = '+str(NbnotMasked), xy=(0.02, 0.8), xycoords='figure fraction', **cp.fig_style.cfont_clck)

    print('     ===> saving figure: '+cfig)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)

    return 0


def ShowBuoysMap( pt, pvlon, pvlat, pvIDs=[], cfig='buoys_RGPS.png', cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1. ):
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

    return _figMap_( pt, pvlon, pvlat, PROJ, cdate=ct, pvIDs=pvIDs, cfig=cfig, ms=ms, ralpha=ralpha, caller='ShowBuoysMap', rzoom=zoom )



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
    #
    cfig = 'debug_interp_'+cname+'_ID'+'%6.6i'%(iID)+'.'+fig_type
    fig = plt.figure(num=1, figsize=(12,5), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes([0.05, 0.13, 0.9, 0.8])
    #plt.axis([ min(vTt)-rdt, max(vTt)+rdt, min(xlat[:,ic])-0.01, max(xlat[:,ic])+0.01])
    #
    func = np.vectorize(dt.utcfromtimestamp)
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
    from math import log
    
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

    plt.savefig(cfig)
    plt.close(1)






def ShowDeformation( pX, pY, pF, cfig='deformation_map.png', cwhat='div', zoom=1,
                     pFmin=-1., pFmax=1., rangeX=None, rangeY=None, unit=None, marker_size=None ):
    '''
    ### Show points, triangle, and quad meshes on the map!
    ###
    ###  * lGeoCoor: True   => we expect degrees for `pX,pY` => geographic (lon,lat) coordinates !
    ###           False  => we expect km or m for `pX,pY` => cartesian coordinates !
    ###
    ###     Specify pX_Q & pY_Q when plotting QuadMesh when IDs are not those of the
    ###     traingle world!
    '''
    from math import log
    
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
    
    #if lGeoCoor:
    #    # Geograhic coordinates (lon,lat)
    #    import cartopy.crs as ccrs
    #    Proj = ccrs.PlateCarree()
    #    vfig = (12*zoom,9*zoom)
    #else:
    # Cartesian coordinates (x,y)
    (xA,xB), (yA,yB), (Lx,Ly), (dx,dy), vfig = _set_fig_axis_( pX, pY, zoom=zoom, rangeX=rangeX, rangeY=rangeY )
    
    fig = plt.figure(num=1, figsize=vfig, facecolor='white')
    
    ddx = dx*Ly/Lx
    if unit:
        ax = plt.axes([1.25*ddx/Lx, 1.4*1.25*dy/Ly, (Lx-2*ddx)/Lx, (Ly-2.2*dy)/Ly], facecolor='0.75')
    else:
        ax = plt.axes([1.25*ddx/Lx,     1.25*dy/Ly, (Lx-2*ddx)/Lx,   (Ly-2*dy)/Ly], facecolor='0.75')
        
    plt.axis([ xA,xB , yA,yB ])

    # Pixel size:
    if not marker_size:
        marker_size = 500.
        if rangeX and rangeY:
            rrm = int( max( abs(rangeX[1]-rangeX[0]) , abs(rangeY[1]-rangeY[0]) ) )
            marker_size = 1710000./rrm
    
    # Showing points:
    #plt.plot( pX, pY, '.', ms=marker_size, color=clPoints, zorder=200) ; #, alpha=0.5)
    plt.scatter( pX, pY, c=pF, s=marker_size, marker='s', cmap=cm, norm=cn )


    if unit:
        # => triggers the colorbar! lilo
        ax2 = plt.axes([0.18,0.03,0.7,0.02])
        #clb = mpl.colorbar.ColorbarBase(ax=ax2, ticks=fa.vc_fld_powlog, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='neither')
        clb = mpl.colorbar.ColorbarBase(ax=ax2, cmap=cm, norm=cn, orientation='horizontal', extend='both')
        clb.set_label(unit, **cfont_clb)
    
    # Adding quadrangles:
    #if len(QuadMesh)>0:
    #    (nbQ,_) = np.shape(QuadMesh)
    #    (nbP,)  = np.shape(pX)
    #    
    #    for jQ in range(nbQ):
    #        vids = QuadMesh[jQ,:] ; # the 4 IDs of the 4 points defining this Quad
    #        if len(pX_Q)>0 and len(pY_Q)>0:
    #            vx, vy = pX_Q[vids], pY_Q[vids]
    #        else:
    #            vx, vy = pX[vids], pY[vids]
    #        vX, vY = np.concatenate([vx[:], vx[0:1]]), np.concatenate([vy[:],vy[0:1]])  ; # need to make an overlap to close the line
    #        plt.plot(vX,vY, col_blu, lw=5*zrat, zorder=100)
    #        plt.fill_between(vX, vY, fc=col_blu, zorder=150, alpha=0.4)
    #        # Indicate quadrangle # in its center:
    #        rmLon, rmLat = np.mean( vx ), np.mean( vy ) ; # Lon,Lat at center of triangle
    #        ax.annotate(str(jQ), (rmLon, rmLat), color='w', fontweight='bold', zorder=160)
    #
    #if len(pnames)>0:
    #    i=0
    #    for cnm in pnames:
    #        ax.annotate(cnm, (pX[i], pY[i]), color=clPNames, fontweight='bold', zorder=500)
    #        i=i+1
    #
    plt.savefig(cfig)
    plt.close(1)

    

def PlotPDFdef( pbinb, pbinc, ppdf, Np=None, name='Divergence', cfig='PDF.png',
                xrng=None, wbin=None, title=None ):
    '''
      * pbinb: vector of the bounds of the bins (x-axis), size = nB+1
      * pbinc: vector of the center of the bins (x-axis), size = nB
      * ppdf:  the PDF, same size as `pbinc`, size = nB

      * wbin:  width of the bin in [day^-1]
    '''
    from math import ceil
    
    nB = len(ppdf)
    if len(pbinc) != nB:
        print('\n *** ERROR ['+caller+'/PlotPDFdef]: wrong size for `pbinc`!'); exit(0)
    if len(pbinb) != nB+1:
        print('\n *** ERROR ['+caller+'/PlotPDFdef]: wrong size for `pbinb`!'); exit(0)

    # Subsampling of labels on the x-axis:
    nx_subsamp = 1
    if xrng and wbin:
        nx_subsamp = int( (xrng[1]-xrng[0])/wbin/10. )
        
    ki = _initStyle_()
        
    width_bin = pbinb[2]-pbinb[1]
        
    fig = plt.figure( num = 1, figsize=(10,9), dpi=None )
    #
    ax1 = plt.axes([0.1, 0.085, 0.85, 0.85])
    #
    #ax1.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    #
    # Y-axis:
    #rinc = 10. ; # => dy of 0.1
    rinc = 100. ; # => dy of 0.01
    ymax = ceil(rinc*np.max(ppdf))/rinc
    ax1.set_ylabel(r'Probability')
    plt.yticks( np.arange(0.,ymax+1./rinc,1./rinc) )
    ax1.set_ylim(0.,ymax)
    #
    # X-axis:
    plt.xlabel(r''+name+' [day$^{-1}$]', color='k')
    #
    if xrng:
        ax1.set_xlim(xrng[0], xrng[1])
        plt.xticks( np.arange(xrng[0], xrng[1]+width_bin, width_bin) )
    else:
        ax1.set_xlim(pbinb[0], pbinb[nB])
        plt.xticks( pbinb[:] )
    #
    if nx_subsamp>1:
        ax_lab = []
        locs, labels = plt.xticks()
        cpt = 0 # with ipct = 1: tick priting will start at y1+dt on x axis rather than y1
        for rr in locs:
            if cpt % nx_subsamp == 0:
                if rr%1.0 == 0.:
                    cr = str(int(rr))  # it's something like 22.0, we want 22 !!!
                else:
                    cr = str(round(rr,6))
                ax_lab.append(cr)
            else:
                ax_lab.append(' ')
            cpt = cpt + 1
        plt.xticks(locs,ax_lab)
        del ax_lab
    
    #plt.plot(pbinc[:], ppdf[:], 'o', color='0.6', zorder=5)
    #
    #plt.bar ( pbinc[:],  ppdf[:], width=width_bin, color='0.6', edgecolor='b', linewidth=2, zorder=10 )
    plt.bar ( pbinc[:],  ppdf[:], width=width_bin, color='0.6', edgecolor=None, linewidth=2, zorder=10 )
    #
    plt.step( pbinb[1:], ppdf[:],  color='k', linewidth=0.5, zorder=10) ; # the enveloppe

    if Np:
        ax1.annotate('N = '+str(Np), xy=(0.72, 0.85), xycoords='figure fraction', **cfont_clock)
    if wbin:
        ax1.annotate('Bin width = '+str(wbin)+r' day$^{-1}$', xy=(0.62, 0.82), xycoords='figure fraction', **cfont_clock)
    if title:
        ax1.annotate(title, xy=(0.5, 0.95), xycoords='figure fraction', ha='center', **cfont_ttl)
        
    plt.savefig(cfig, dpi=100, orientation='portrait', transparent=False)
    plt.close(1)

    return 0





def PlotMesh( pcoor_trg, Ys, Xs, isrc_msh, vnames=['P1','P2','P3','P4'], fig_name='mesh.png', pcoor_extra=(-999.,-999.), label_extra=None ):
    '''
    isrc_msh: 2D integer array of shape (4,2)
    wghts:    1D real array of shape (4,)
    pcoor_extra: just an extra point to show on the figure....
    '''
    (yT,xT)                             = pcoor_trg
    [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = isrc_msh[:,:]
    #
    fig = plt.figure(num = 1, figsize=[14,10], facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.09, 0.07, 0.6, 0.9])

    plt.plot( [     xT    ], [     yT    ], marker='o', ms=15, color='k', label='target point' ) ; # target point !
        
    for i in range(4):
        [ja,ia] = isrc_msh[i,:]
        [jb,ib] = isrc_msh[(i+1)%4,:]
        # Lines to see the mesh:
        plt.plot( [ Xs[ja,ia],Xs[jb,ib] ], [ Ys[ja,ia],Ys[jb,ib] ], linestyle='-', color='k', marker=None, ms=10, label=None)
        # The point:
        plt.plot( [ Xs[ja,ia] ], [ Ys[ja,ia] ], marker='o', ms=10, label=vnames[i] )

    if pcoor_extra!=(-999.,-999.):
        (yE,xE) = pcoor_extra
        plt.plot( [     xE    ], [     yE    ], marker='+', ms=20, color='0.5', label=label_extra ) ; # target point !

    ax1.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), fancybox=True)
    plt.savefig(fig_name, dpi=100, transparent=False)
    plt.close(1)
