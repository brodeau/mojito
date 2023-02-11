
from sys import argv, exit
from os import path ; #, environ, mkdir
import numpy as np

#from re import split
from netCDF4 import Dataset

#from climporn import epoch2clock, clock2epoch
#import mojito as mjt

tunits_default = 'seconds since 1970-01-01 00:00:00'


def chck4f( cf ):
    if not path.exists(cf):
        print(' ERROR [chck4f()]: file '+cf+' does not exist!')
        exit(0)


def InspectLoadData( cfile, plistVar, iverbose=0 ):
    '''
       Open & inspect the RGPS NetCDF input file and load the raw data
    '''
    #
    chck4f(cfile)
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        if not 'points' in list_dim:
            print(' ERROR [InspectLoadData()]: no dimensions `points` found into input file!'); exit(0)

        list_var = list( id_in.variables.keys() )
        if iverbose>0: print(' ==> variables: ', list_var, '\n')
        #
        for cv in plistVar:
            if not cv in list_var:
                print(' ERROR [InspectLoadData()]: no variable `'+cv+'` found into input file!'); exit(0)

        nP = id_in.dimensions['points'].size
        if iverbose>0: print('\n *** Total number of points in the file = ', nP)

        # Time records:
        ctunits = id_in.variables['time'].units
        if not ctunits == tunits_default:
            print(" ERROR [InspectLoadData()]: we expect '"+tunits_default+
                  "' as units for the time record vector, yet we have: "+ctunits)
            exit(0)
        ztime = id_in.variables['time'][:]

        # Coordinates:
        zx    = id_in.variables['x'][:]
        zy    = id_in.variables['y'][:]
        zlon  = id_in.variables['lon'][:]
        zlat  = id_in.variables['lat'][:]

        # Buoys' IDs:
        kBIDs    = np.zeros(nP, dtype=int)
        kBIDs[:] = id_in.variables['index'][:]

    zlon[:] = np.mod(zlon, 360.) ; # Longitudes in the [0:360] frame...

    return nP, ztime, zx, zy, zlon, zlat, kBIDs



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

