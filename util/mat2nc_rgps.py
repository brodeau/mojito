#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

import datetime as dt
import glob
from os import path

from cartopy.crs import PlateCarree, NorthPolarStereo
import numpy as np
from scipy.io import loadmat

import mojito as mjt

srs_longlat = PlateCarree()
srs_rgps = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)

#idir = '/data/gcm_setup/data/RGPS_Kwok_98/Ramplab'
idir = './mat'
#
#dst_dir = '/data/gcm_setup/data/RGPS_Kwok_98'
dst_dir = '.'



if __name__ == '__main__':

    ifiles = sorted(glob.glob(f'{idir}/RGPS_*_traj.mat'))
    
    for ifile in ifiles:
        
        dst_file = f'{dst_dir}/{path.basename(ifile).replace(".mat", ".nc4")}'
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


        kk = mjt.SaveRGPStoNC( dst_file, time_all, buoy_ids, y_all, x_all, lat_all, lon_all, pqual=qual_all, Nstrm=NbStreams, pSid=strm_all )
        #kk = mjt.SaveRGPStoNC( dst_file, time_all, buoy_ids, y_all, x_all, lat_all, lon_all, qual_all )

