
from sys import argv, exit
from os import path ; #, environ, mkdir
import numpy as np

#from re import split
from netCDF4 import Dataset

#from climporn import epoch2clock, clock2epoch
#import mojito as mjt

tunits_default = 'seconds since 1970-01-01 00:00:00'


def ncSaveCloudBoys( cf_out, ptime, pIDs, pY, pX, pLat, pLon, tunits=tunits_default ):
    '''
    '''
    print('\n *** [ncSaveCloudBoys]: About to generate file: '+cf_out+' ...')
    (Nt,) = np.shape(ptime)
    (Nb,) = np.shape(pIDs)
    if np.shape(pY)!=(Nt,Nb) or np.shape(pX)!=(Nt,Nb) or np.shape(pLat)!=(Nt,Nb) or np.shape(pLon)!=(Nt,Nb):
        print('ERROR [ncSaveCloudBoys]: one of the 2D arrays has a wrong shape!!!')
        exit(0)
    
    f_out = Dataset(cf_out, 'w', format='NETCDF4')

    # Dimensions:
    cd_time = 'time'
    cd_buoy = 'buoy'
    f_out.createDimension(cd_time, None)
    f_out.createDimension(cd_buoy, Nb  )

    # Variables:
    v_time  = f_out.createVariable(cd_time,     'i4',(cd_time,))
    v_buoy  = f_out.createVariable(cd_buoy,     'i4',(cd_buoy,))
    v_bid   = f_out.createVariable('id_buoy',   'i4',(cd_buoy,))
    v_latb  = f_out.createVariable('latitude' , 'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)
    v_lonb  = f_out.createVariable('longitude', 'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)
    v_y     = f_out.createVariable('y_pos' ,    'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)
    v_x     = f_out.createVariable('x_pos',     'f4',(cd_time,cd_buoy,), zlib=True, complevel=9)

    v_time.units = tunits
    v_bid.units  = 'ID of buoy'
    v_latb.units = 'degrees north'
    v_lonb.units = 'degrees south'
    v_y.units    = 'km'
    v_x.units    = 'km'

    v_buoy[:] = np.arange(Nb,dtype='i4')
    v_bid[:]  = pIDs[:]

    for jt in range(Nt):
        v_time[jt]   = ptime[jt]
        v_latb[jt,:] =  pLat[jt,:]
        v_lonb[jt,:] =  pLon[jt,:]
        v_y[jt,:]    =    pY[jt,:]
        v_x[jt,:]    =    pX[jt,:]

    f_out.About  = 'RGPS sea-ice drift data (Kwok, 1998)'
    f_out.Author = 'Generated with `'+path.basename(argv[0])+'` of `mojito` (L. Brodeau, 2022)'
    f_out.close()
    print('      ===> '+cf_out+' saved!')

    return 0
