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
import lbrgps   as lbr

l_control_IDs = False

list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

idebug = 1

narg = len(argv)
if not narg in [3]:
    print('Usage: '+argv[0]+' <file_RGPS.nc> <buoy_ID_to_follow>')
    exit(0)
cf_in  = argv[1]
id_fl  = int(argv[2])



# Getting time info and time step from input model file:
vv = split('_', path.basename(cf_in))
cdt2 = vv[-2]
cdt1 = vv[-3]
print('\n *** Date #1 = '+cdt1)
print('\n *** Date #2 = '+cdt2)
#vv = split('-|_', path.basename(cf_in))
#cyear = vv[-4]
#print('\n *** Year = '+cyear)
#exit(0)
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
            exit(0)


vIDs = nmp.arange(id_min, id_max+1) ; # this should be all the buoys IDs
Nb   = len(vIDs)

print("\n *** There are "+str(Nb)+" buoys to follow: ID "+str(id_min)+" to ID "+str(id_max)+" !")

if not id_fl in vIDs:
    print("\n ERROR: buoy ID #'+str(id_fl)+' does not exist in the file!")

print("\n *** We are going to follow buoy with ID: "+str(id_fl))

# Indexes of "points" that deal with this buoy:
idx_buoy, = nmp.where( vindex == id_fl )
Nt = len(idx_buoy)


# What is going to be plotted:

vIDs = [ id_fl ]
Nbf = len(vIDs)

vtimb = nmp.zeros( Nt   , dtype=int)
vlonb = nmp.zeros((Nt,Nbf))
vlatb = nmp.zeros((Nt,Nbf))

vtimb[:]      = vtime[ idx_buoy ]
if Nbf>1:
    vlonb[:,:Nbf] =  vlon[ idx_buoy ]
    vlatb[:,:Nbf] =  vlat[ idx_buoy ]
else:
    vlonb[:,0] =  vlon[ idx_buoy ]
    vlatb[:,0] =  vlat[ idx_buoy ]



ctime = nmp.zeros(Nt, dtype='U32')

print("     ===> it has "+str(Nt)+" time records")
for jt in range(Nt):
    ctime[jt] = cp.epoch2clock(vtimb[jt])
    if l_control_IDs: print("       jt="+str(jt)+" => time = "+ctime[jt])
    

# Figures:
    
cFigPref = 'buoy_'+'%6.6i'%(id_fl)

kf = lbr.ShowBuoysMap_Trec( vtimb, vlonb, vlatb, pvIDs=vIDs, cnmfig=cFigPref, ms=20, ralpha=1 )
