#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

import datetime as dt

import numpy as np
#from scipy.io import loadmat
from netCDF4 import Dataset

import pandas as pd

from cartopy.crs import NorthPolarStereo

list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

srs_nextsim = NorthPolarStereo(central_longitude=-45, true_scale_latitude=60)
srs_rgps    = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)


#ifile = 'RGPS_2006-12-03_2007-06-01_traj.mat'
ifile = '/MEDIA/data/data/RGPS_Kwok_98/RGPS_2006-12-03_2007-06-01_traj.nc4'
n_streams = 15

# load file
# netCDF rather than matlab....
id_in = Dataset(ifile)
#m     = loadmat(ifile)

list_dim = list( id_in.dimensions.keys() )
if not 'points' in list_dim:
    print(' ERROR: no dimensions `points` found into input file!'); exit(0)
#
list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
for cv in list_expected_var:
    if not cv in list_var:
        print(' ERROR: no variable `'+cv+'` found into input file!'); exit(0)
#
Np0 = id_in.dimensions['points'].size
print(' *** Number of provided virtual buoys = ', Np0,'\n')


streams = []
#for i in range(len(m['out'][0])):
for i in range(Np0):
    print(i)
    stream = m['out'][0][i]

    x    = [j[0] for j in stream['trajectories'][0][0]['x_map' ][0]]
    y    = [j[0] for j in stream['trajectories'][0][0]['y_map' ][0]]
    time = [j[0] for j in stream['trajectories'][0][0]['time'  ][0]]
    year = [j[0] for j in stream['trajectories'][0][0]['year'  ][0]]
    q    = [j[0] for j in stream['trajectories'][0][0]['q_flag'][0]]
    
    d    = [np.array([dt.datetime(yi,1,1) + dt.timedelta(ti) for (yi,ti) in zip(y_vec, t_vec)])
                    for (y_vec,t_vec) in zip(year, time)]

    x_new = []
    y_new = []
    # convert to meters in neXtSIM projection
    for x0,y0 in zip(x,y):
        x1,y1,_ = srs_nextsim.transform_points(srs_rgps, x0, y0).T * 1000.
        x_new.append(x1)
        y_new.append(y1)

    streams.append({
    'x': x_new,
    'y': y_new,
    'd': d,
    'q': q,
    })

pairs = []
for stream in streams:
    for x_traj, y_traj, d_traj in zip(stream['x'], stream['y'], stream['d']):
        for x1, y1, d1, x2, y2, d2 in zip(
                x_traj[:-1], y_traj[:-1], d_traj[:-1], x_traj[1:], y_traj[1:], d_traj[1:]):
            pairs.append([x1, x2, y1, y2, d1, d2])   
df = pd.DataFrame(pairs, columns=['rx0', 'rx1', 'ry0', 'ry1', 'rd0', 'rd1'])
df.to_pickle(ifile.replace('.mat', '_df.npz'))
