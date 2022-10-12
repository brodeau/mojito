#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import exit, argv
from os import path
import numpy as nmp
import pandas as pd

from re import split

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid

from datetime import datetime as dt

import climporn as cp
import lbrgps   as lbr

list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

idebug = 1

rtol_sec = 3600.*12. # tolerance in seconds, to consider it's the "same time" !!! :(


narg = len(argv)
if not narg in [3]:
    print('Usage: '+argv[0]+' <file_RGPS.nc> <buoy_ID_to_follow (single: `ID`, or range: `ID1,ID2`>')
    exit(0)
cf_in = argv[1]
cf_bf = argv[2]

l_id_range = (',' in cf_bf)
if l_id_range:
    vv = split(',', path.basename(cf_bf))
    ib1 = int(vv[0])
    ib2 = int(vv[1])
    if ib2<=ib1: print("ERROR: your buoy ID range makes no sense!"); exit(0)
    vid_fl  = nmp.arange(ib1,ib2+1)
else:
    vid_fl  = [ int(argv[2]) ]
    
#print("\n *** List of buoys' IDs to follow: ", vid_fl[:])

Nbf = len(vid_fl)

print("\n *** There are "+str(Nbf)+" buoys to follow: ID "+str(ib1)+" to ID "+str(ib2)+" !")
#print( vid_fl )

# Getting time info and time step from input model file:
#vv = split('_', path.basename(cf_in))
#cdt2 = vv[-2]
#cdt1 = vv[-3]
#print('\n *** Date #1 = '+cdt1)
#print('\n *** Date #2 = '+cdt2)



# Opening and inspecting input file
cp.chck4f(cf_in)
#
id_in    = Dataset(cf_in)
#
list_dim = list( id_in.dimensions.keys() ) ; #print(' ==> dimensions: ', list_dim, '\n')
if not 'points' in list_dim:
    print(' ERROR: no dimensions `points` found into input file!'); exit(0)
#
list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
for cv in list_expected_var:
    if not cv in list_var:
        print(' ERROR: no variable `'+cv+'` found into input file!'); exit(0)
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






# Of all the buoys selected in range we are going to scan them all to see which one has the longest record!


ib_longest = -1 ; Ntmax = 0 ; idx_buoy_longest = [] ; ID_longest = -1 ; vt_ref = []
ib = 0
for id_fl in vid_fl:
    if not id_fl in vid_fl:  print("\n ERROR: buoy ID #'+str(id_fl)+' does not exist in the file!"); exit(0)        
    idx_buoy, = nmp.where( vindex == id_fl ) ; # Indexes of "points" that deal with this buoy:
    Nt = len(idx_buoy)
    if Nt > Ntmax:
        Ntmax      = Nt
        ID_longest = id_fl
        ib_longest = ib
        idx_buoy_longest = idx_buoy
        vt_ref = vtime[idx_buoy]
    ib = ib + 1

print("\n *** Buoy with the longest record has ID #"+str(ID_longest)+" !")
print("     ==> number of records: Nmax = "+str(Ntmax))


# Now, the real shit:
vmskID = nmp.zeros( Nbf,        dtype=int)
IX     = nmp.zeros((Ntmax,Nbf), dtype=int) - 999 ; #fixme!  2000 is max number of records for a given buoy!
IM     = nmp.zeros((Ntmax,Nbf), dtype=int)       ; #fixme!  2000 is max number of records for a given buoy!
vib    = [] ; # will store the valid `ib` indices...
ib = 0
for id_fl in vid_fl:
    idx_buoy, = nmp.where( vindex == id_fl )
    Nt = len(idx_buoy)
    # We must ensure that the time records are pretty much the same as `vt_ref` !!!
    vt = vtime[idx_buoy]    
    if abs(vt_ref[0]-vt[0]) < rtol_sec and abs(vt_ref[Nt-1]-vt[Nt-1]) < rtol_sec:    
        if idebug>1: print("   * buoy #"+str(id_fl)+" works!")
        vmskID[ib] = 1
        IX[0:Nt,ib] = idx_buoy
        IM[0:Nt,ib] = 1
        vib.append(ib)
    ib = ib + 1

IX   = nmp.ma.masked_where( IM==0,     IX )
vIDs = nmp.ma.MaskedArray.compressed( nmp.ma.masked_where( vmskID==0, vid_fl ) ) ; # valid remaining IDs...

Nbv = len(vIDs) ; # update Nbv !
if Nbv != len(vib): print('ERROR ZZ!'); sys.exit(0)
if Nbv < 2:         print("ERROR: not enough valid buoys!"); exit(0)
print("\n *** There are now "+str(Nbv)+" VALID buoys to follow!")

# What is going to be plotted:
vtimb = nmp.zeros( Ntmax   , dtype=int)
vlonb = nmp.zeros((Ntmax,Nbv))
vlatb = nmp.zeros((Ntmax,Nbv))
vmskb = nmp.zeros((Ntmax,Nbv) , dtype=int)

vtimb[:]      = vtime[ idx_buoy_longest ]

for iv in range(Nbv):
    ib = vib[iv]
    print("ib = ", ib)
    idx_buoy = nmp.ma.MaskedArray.compressed(IX[0:Nt,ib])
    print("ib = ", ib, " => idx_buoy =", idx_buoy)
    nl = len(idx_buoy)
    vlonb[:nl,iv] =  vlon[ idx_buoy ]
    vlatb[:nl,iv] =  vlat[ idx_buoy ]
    vmskb[:nl,iv] = 1

del vlon, vlat, vtime

vlonb = nmp.ma.masked_where( vmskb==0, vlonb )
vlatb = nmp.ma.masked_where( vmskb==0, vlatb )


# Figures:
    
#cFigPref = 'buoys_'+'%6.6i'%(vIDs)
cFigPref = 'buoys'

kf = lbr.ShowBuoysMap_Trec( vtimb, vlonb, vlatb, pvIDs=vIDs, cnmfig=cFigPref )


