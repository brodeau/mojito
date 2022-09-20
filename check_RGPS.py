#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

import sys
from os import path
import numpy as nmp
import pandas as pd

from re import split

from netCDF4 import Dataset

#import csv

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid

#from calendar import isleap
from datetime import datetime as dt

import climporn as cp

l_control_IDs = False

list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

idebug = 1

color_top = 'w'
clr_yellow = '#ffed00'

# Projection:
#vp =  ['Arctic', 'stere', -60., 40., 122., 57.,    75.,  -12., 10., 'h' ]  # Big Arctic + Northern Atlantic
vp =  ['Arctic', 'stere', -80., 68., 138.5, 62.,    90.,  -12., 10., 'h' ]  # North Pole Arctic (zoom)

rfig_fact = 1.
vfig_size = [ 7.54*rfig_fact, 7.2*rfig_fact ]

vsporg = [0., 0., 1., 1.]
col_bg = '#041a4d'

# Defaults:
#lknown = True
#rfact_zoom = 1.
#l_show_cb = False
#l_show_nm = False
#l_show_clock=False ; x_clock=0. ; y_clock=0.
#l_scientific_mode = False
#l_show_ttl = False
#vcb = [0.15, 0.96, 0.8, 0.02]
font_rat = 1.

#l_show_msh = False

pt_sz_track = 1

fig_type='png'
rDPI = 150

#vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
#vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]


narg = len(sys.argv)
if not narg in [2]:
    print('Usage: '+sys.argv[0]+' <file_RGPS.nc>')
    sys.exit(0)

cf_in  = sys.argv[1]
#cf_in = sys.argv[2]
#vv = split(',',sys.argv[3])
#cf_in = vv[0]
#cv_mod = vv[1]
#cnfig  = sys.argv[4]
#
#cf_lsm = sys.argv[5]
#
# Subsampling in time...
#itsubs = 1
#if narg == 7 :
#    itsubs = int(sys.argv[6])



# Getting time info and time step from input model file:
vv = split('_', path.basename(cf_in))
cdt2 = vv[-2]
cdt1 = vv[-3]
print('\n *** Date #1 = '+cdt1)
print('\n *** Date #2 = '+cdt2)
#vv = split('-|_', path.basename(cf_in))
#cyear = vv[-4]
#print('\n *** Year = '+cyear)
#sys.exit(0)
#vm = vmn
#if isleap(yr0): vm = vml


#dir_conf = path.dirname(cf_in)
#if dir_conf == '':  dir_conf = '.'
#print('\n *** dir_conf =',dir_conf,'\n')



# Opening and inspecting input file
cp.chck4f(cf_in)

id_in    = Dataset(cf_in)
#
list_dim = list( id_in.dimensions.keys() ) ; #print(' ==> dimensions: ', list_dim, '\n')
if not 'points' in list_dim:
    print(' ERROR: no dimensions `points` found into input file!'); sys.exit(0)
#
list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
for cv in list_expected_var:
    if not cv in list_var:
        print(' ERROR: no variable `'+cv+'` found into input file!'); sys.exit(0)
#
Np = id_in.dimensions['points'].size
print(' *** Number of provided virtual buoys = ', Np)


vlon  = id_in.variables['lon'][:]
vlat  = id_in.variables['lat'][:]

rlat_min = nmp.min(vlat)
print('     ==> Southernmost latitude = ', rlat_min)

rlat_max = nmp.max(vlat)
print('     ==> Northernmost latitude = ', rlat_max)

#for jj in range( 0,Np,10 ):  print('  ', jj, vlat[jj])
#for ji in range( 0,Np,10 ):  print('  ', ji, vlon[ji])

vtime = id_in.variables['time'][:]

vindex = nmp.zeros(Np, dtype=int)
vindex[:] = id_in.variables['index'][:]

id_in.close()


## How many buoys (IDs) are there?
id_min = nmp.min(vindex)
id_max = nmp.max(vindex)
#print(" Min, max ID:", id_min, id_max )

if l_control_IDs:
    for jid in range(id_min, id_max+1):
        if jid in vindex[:]:
            vIDs.append(jid)
        else:
            print(' there is no buoy with ID:', jid, '!')
            sys.exit(0)


vIDs = nmp.arange(id_min, id_max+1) ; # this should be all the buoys IDs
Nb   = len(vIDs)

print("\n *** There are "+str(Nb)+" buoys to follow: ID "+str(id_min)+" to ID "+str(id_max)+" !\n")


# Now we must id all different times that there is:
#vdate = pd.unique(vtime)

vdates,idx_d  = nmp.unique( vtime, return_index=True)

vtmp = nmp.sort(vdates)
if not ( nmp.sum(nmp.abs(vdates-vtmp)) == 0.) :
    print("Date vector is not sorted...")
    sys.exit(0)
del vtmp

Nt = len(vdates)
print(" *** There are "+str(Nt)+" different time snapshots !")
print("     ==> going from '"+str(dt.fromtimestamp(vdates[0]))+" to "+str(dt.fromtimestamp(vdates[Nt-1]))+"\n")


# Construct vector of daily dates (format: '2007-12-02'):
cdates = nmp.zeros(Nt, dtype='U10')
for jt in range(Nt):
    cdates[jt] = dt.fromtimestamp(vdates[jt]).strftime("%Y-%m-%d")


# Keep unique dates (daily then):
CDD = nmp.unique( cdates )
Nd  = len(CDD)     
print("\n *** There are "+str(Nd)+" different days !")
print("     ==> going from '"+CDD[0]+" to "+CDD[Nd-1]+"\n")
#sys.exit(0)


print("\n *** Building a copy of time axis as strings of daily dates ...")
ctdays = nmp.zeros(Np, dtype='U10')  # for all the points we need to know the day => Np !!!
for jp in range(Np):
    ctdays[jp] = dt.fromtimestamp(vtime[jp]).strftime("%Y-%m-%d")
    

print(ctdays[:])
sys.exit(0)
    

#### FIGURE ####

#id_floats2follow = nmp.arange(0,Np,10)
#print(id_floats2follow)

PROJ = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
               resolution=vp[9], area_thresh=1000., projection='stere', \
               lat_0=vp[6], lon_0=vp[7], epsg=None)



#for jt in range(Nt):
for jd in range(Nd):

    ctd = CDD[jd]

    print('\n *** Plotting for day = '+ctd)
    
    #cfig = 'test_'+'%4.4i'%(jt+1)+'.'+fig_type
    cfig = 'test_'+ctd+'.'+fig_type

    #sys.exit(0)
    

    # Daily current date:
    #ctd = dt.fromtimestamp(vdates[jt]).strftime("%Y-%m-%d")

    #print('\n *** Plotting for time ='+str(dt.fromtimestamp(vdates[jt])))
    #print('     => widening to day: '+ctd)

    # List of snapshots that are in this current day:
    idx_t, = nmp.where( cdates[:] == ctd )

    vsnaps = vtime[idx_t]
    
    print(vsnaps)
    sys.exit(0)

    
    

    
    fig = plt.figure(num=1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
    ax  = plt.axes(vsporg, facecolor = col_bg)

    x0,y0 = PROJ(vlon[idx_t],vlat[idx_t])
    #ct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=pt_sz_track )
    ct = plt.scatter(x0, y0, color='w',  marker='.', s=pt_sz_track )

    PROJ.drawcoastlines(linewidth=0.5)
    PROJ.fillcontinents(color='grey') #, alpha=0)
    #PROJ.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    #PROJ.drawmapboundary()
    PROJ.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1], linewidth=0.3)
    PROJ.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0], linewidth=0.3)

    print(' Saving figure: '+cfig)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)



    
#for jp in range(0,Np,10):
#    
#    #print('    jp =', jp)
#    #ct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=pt_sz_track )
#    #ct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=pt_sz_track )
#    #ct = plt.plot(x0, y0) ; #, marker='.', s=pt_sz_track )
#
#    if vIDs[jp] in id_floats2follow:
#        x0,y0 = PROJ(vlon[jp],vlat[jp])
#        #ct = plt.scatter(x0, y0, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=pt_sz_track )
#        ct = plt.scatter(x0, y0, color='w',  marker='.', s=pt_sz_track )
        
#x0,y0 = PROJ(Xlon,Xlat)

#f = PROJ.pcolor(x0, y0, XIFR, cmap = pal_ice, norm = nrm_ice) #, interpolation='none')

#if lshow_ice: XFLD = nmp.ma.masked_where(XIFR >= 0.2, XFLD)
#XFLD = nmp.ma.masked_where(XMSK <= 0.1, XFLD)

#if lshow_ice:   XIFR = nmp.ma.masked_where(XIFR  < r_oi_thr, XIFR)

#if l_only_over_ice:   XFLD = nmp.ma.masked_where(XIFR <  r_oi_thr, XFLD)


#ft = PROJ.pcolormesh(x0, y0, XFLD, cmap = pal_fld, norm = nrm_fld, shading="gouraud" )
#ft = PROJ.pcolormesh(x0, y0, XFLD, cmap = pal_fld, norm = nrm_fld )
#ft = PROJ.pcolor(x0, y0, XFLD, cmap = pal_fld, norm = nrm_fld )



# ----------- Color bar for field -----------
#ax2 = plt.axes([0.64, 0.965, 0.344, 0.018])
#clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=nrm_fld, orientation='horizontal', extend='both')
#cb_labs = []
#cpt = 0
#for rr in vc_fld:
#    if cpt % cb_jump == 0 or ( (tmin == -tmax) and (int(rr) == 0 ) ):
#        if df >= 1.: cb_labs.append(str(int(rr)))
#        if df <  1.: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))
#    else:
#        cb_labs.append(' ')
#    cpt = cpt + 1
#clb.ax.set_xticklabels(cb_labs, **cfont_clb)
#clb.set_label(cunit, **cfont_clb)
#clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color
#clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor
#clb.ax.tick_params(which = 'minor', length = 2, color = color_top_cb )
#clb.ax.tick_params(which = 'major', length = 4, color = color_top_cb )
#del clb

#ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.3, 0.01), xycoords='figure fraction', **cfont_clock)

#ry0 = 0.78
#ax.annotate(cxtra_info1, xy=(0.02, ry0+0.05), xycoords='figure fraction', **cfont_titl1)

#if ctitle != "":
#    ax.annotate(ctitle, xy=(0.03, ry0-0.01), xycoords='figure fraction', **cfont_titl2)



sys.exit(0)





params = { 'font.family':'Open Sans',
           'font.weight':    'normal',
           'font.size':       int(12.*font_rat),
           'legend.fontsize': int(22.*font_rat),
           'xtick.labelsize': int(18.*font_rat),
           'ytick.labelsize': int(18.*font_rat),
           'axes.labelsize':  int(15.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
cfont_date = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(12.*font_rat), 'color':'w' }
cfont_mail = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_cnfn = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(35.*font_rat), 'color':'w' }
cfont_axis  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
cfont_ttl = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(25.*font_rat), 'color':color_top }
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(18.*font_rat), 'color':color_top }

# Colormaps for fields:
pal_fld = cp.chose_colmap(cpal_fld)
if l_log_field:
    norm_fld = colors.LogNorm(  vmin = tmin, vmax = tmax, clip = False)
if l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm = cp.chose_colmap('land_dark')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

pal_filled = cp.chose_colmap('gray_r')
norm_filled = colors.Normalize(vmin = 0., vmax = 0.1, clip = False)

fsize  = ( rh, rh*yx_ratio )
vc_fld = nmp.arange(tmin, tmax + df, df)













#cp.chck4f(cf_lsm)
#cnmsk = 'tmask'
#print('\n *** Reading "'+cnmsk+'" in meshmask file...')
#id_lsm = Dataset(cf_lsm)
#nb_dim = len(id_lsm.variables[cnmsk].dimensions)
#Ni = id_lsm.dimensions['x'].size
#Nj = id_lsm.dimensions['y'].size
#if i2 == 0: i2 = Ni
#if j2 == 0: j2 = Nj
#if nb_dim == 4: XMSK  = id_lsm.variables[cnmsk][0,0,j1:j2,i1:i2] ; # t, y, x
#if nb_dim == 3: XMSK  = id_lsm.variables[cnmsk][0,  j1:j2,i1:i2] ; # t, y, x
#if nb_dim == 2: XMSK  = id_lsm.variables[cnmsk][    j1:j2,i1:i2] ; # t, y, x
#if l_show_msh:
#    Xlon = id_lsm.variables['glamu'][0,j1:j2,i1:i2]
#    Xlat = id_lsm.variables['gphiv'][0,j1:j2,i1:i2]
#id_lsm.close()
#print('      done.')


