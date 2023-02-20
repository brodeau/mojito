#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

import datetime as dt
import glob
import os

from cartopy.crs import PlateCarree, NorthPolarStereo
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy.io import loadmat

srs_longlat = PlateCarree()
srs_rgps = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)

idir = '/data/gcm_setup/data/RGPS_Kwok_98/Ramplab'
dst_dir = '/data/gcm_setup/data/RGPS_Kwok_98'
#dst_dir = '.'

ifiles = sorted(glob.glob(f'{idir}/RGPS_*_traj.mat'))



for ifile in ifiles:
    
    dst_file = f'{dst_dir}/{os.path.basename(ifile).replace(".mat", ".nc4")}'
    print(ifile, dst_file)
    
    # load file
    m = loadmat(ifile)

    date0 = dt.datetime(1970,1,1)

    x_all, y_all, time_all, qual_all, buoy_ids = [], [], [], [], []
    strm_all = []
    buoy_id = 0


    data0 = m['out']

    print(' Shape of `data0` =',np.shape(data0))
    (_,NbStreams) = np.shape(data0)

    print( '   => NbStreams =',NbStreams)

    print(len(m['out'][0]))
    

    for jS in range(NbStreams):

        id_stream = jS + 1

        print(' *** Doing stream #',id_stream)

        stream = m['out'][0][jS]
        
        x    = [ j[0] for j in stream['trajectories'][0][0]['x_map'][0]  ]
        y    = [ j[0] for j in stream['trajectories'][0][0]['y_map'][0]  ]
        time = [ j[0] for j in stream['trajectories'][0][0]['time'][0]   ]
        year = [ j[0] for j in stream['trajectories'][0][0]['year'][0]   ]
        q    = [ j[0] for j in stream['trajectories'][0][0]['q_flag'][0] ]
        d    = [  np.array([(dt.datetime(yi,1,1) + dt.timedelta(ti) - date0).total_seconds() for (yi,ti) in zip(y_vec, t_vec)])
                 for (y_vec,t_vec) in zip(year, time) ]
        
        x_all.append(np.hstack(x))
        y_all.append(np.hstack(y))
        time_all.extend(d)
        qual_all.extend(q)
        
        for j in stream['trajectories'][0][0]['x_map'][0]:
            buoy_ids.append( np.ones(j[0].shape) * buoy_id   )
            strm_all.append( np.zeros(j[0].shape) + id_stream )
            buoy_id += 1

    ### for jS in range(NbStreams)

            
    x_all, y_all, time_all, qual_all, buoy_ids, strm_all = [ np.hstack(i) for i in [ x_all, y_all, time_all, qual_all, buoy_ids, strm_all ] ]
    lon_all, lat_all, _ = srs_longlat.transform_points(srs_rgps, x_all*1000, y_all*1000).T    

    
    with Dataset(dst_file, 'w') as ds:
        
        points  = ds.createDimension( "points", lon_all.size )
        streams = ds.createDimension( "streams", NbStreams   )
        
        time   = ds.createVariable("time", "f8", ("points",), zlib=True, complevel=9)
        x      = ds.createVariable("x", "f8", ("points",), zlib=True, complevel=9)
        y      = ds.createVariable("y", "f8", ("points",), zlib=True, complevel=9)
        lon    = ds.createVariable("lon", "f8", ("points",), zlib=True, complevel=9)
        lat    = ds.createVariable("lat", "f8", ("points",), zlib=True, complevel=9)
        qual   = ds.createVariable("q_flag", "byte", ("points",), zlib=True, complevel=9)
        index  = ds.createVariable("index", "long", ("points",), zlib=True, complevel=9)
        idstream = ds.createVariable("idstream", "byte", ("points",), zlib=True, complevel=9)
        streams = ds.createVariable("streams", "byte", ("streams",) )

        time[:] = time_all
        time.setncattr('units', "seconds since 1970-01-01 00:00:00")
        time.setncattr('standard_name', "time")
        time.setncattr('long_name', "time of virtual buoy")
        time.setncattr('calendar', "standard")

        x[:] = x_all
        x.setncattr('units', "km")
        x.setncattr('long_name', "X-coordinate in Polar Stereographic projection, lon_0=-45, lat_ts=70")
        x.setncattr('standard_name', "x_coordinate")

        y[:] = y_all
        y.setncattr('units', "km")
        y.setncattr('long_name', "Y-coordinate in Polar Stereographic projection, lon_0=-45, lat_ts=70")
        y.setncattr('standard_name', "y_coordinate")

        lon[:] = lon_all
        lon.setncattr('units', "degrees_easth")
        lon.setncattr('long_name', "longitude coordinate of virtual buoy")
        lon.setncattr('standard_name', "longitude")

        lat[:] = lat_all
        lat.setncattr('units', "degrees_north")
        lat.setncattr('long_name', "latitude coordinate of virtual buoy")
        lat.setncattr('standard_name', "latitude")

        qual[:] = qual_all
        qual.setncattr('long_name', "Quality flag of virtual buoy")

        index[:] = buoy_ids
        index.setncattr('long_name', "Id of virtual buoy")

        idstream[:] = strm_all
        idstream.setncattr('long_name', "ID of the stream the virtual buoy belongs to")        

        streams[:] = np.arange(NbStreams) + 1
        streams.setncattr('long_name', "IDs of existing streams")
        
        ds.setncattr('title', 'RGPS trajectories')
        ds.setncattr('reference', 'Kwok, Ronald. “The RADARSAT Geophysical Processor System.” (1998).')    

    print(' *** File '+dst_file+' generated!\n')
