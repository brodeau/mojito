
from sys import argv, exit
from os import path ; #, environ, mkdir
import numpy as np

#from re import split
from netCDF4 import Dataset

#from climporn import epoch2clock, clock2epoch

from .util import chck4f, ConvertGeo2CartesianNPSkm

tunits_default = 'seconds since 1970-01-01 00:00:00'

def LoadDataRGPS( cfile, plistVar, iverbose=0 ):
    '''
       Open & inspect the RGPS NetCDF input file and load the raw data
    '''
    #
    chck4f(cfile)
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        if not 'points' in list_dim:
            print(' ERROR [LoadDataRGPS()]: no dimensions `points` found into input file!'); exit(0)

        list_var = list( id_in.variables.keys() )
        if iverbose>0: print(' ==> variables: ', list_var, '\n')
        #
        for cv in plistVar:
            if not cv in list_var:
                print(' ERROR [LoadDataRGPS()]: no variable `'+cv+'` found into input file!'); exit(0)

        nP = id_in.dimensions['points'].size
        if iverbose>0: print('\n *** Total number of points in the file = ', nP)

        # Time records:
        ctunits = id_in.variables['time'].units
        if not ctunits == tunits_default:
            print(" ERROR [LoadDataRGPS()]: we expect '"+tunits_default+
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

    return nP, ztime, zy, zx, zlat, zlon, kBIDs



def GetModelGrid( fNCmeshmask ):

    chck4f( fNCmeshmask)
    
    # Reading mesh metrics into mesh-mask file:
    with Dataset(fNCmeshmask) as id_mm:
        kmaskt = id_mm.variables['tmask'][0,0,:,:]
        zlonF  = id_mm.variables['glamf'][0,:,:]
        zlatF  = id_mm.variables['gphif'][0,:,:]
        zlonT  = id_mm.variables['glamt'][0,:,:]
        zlatT  = id_mm.variables['gphit'][0,:,:]
        ze1T   = id_mm.variables['e1t'][0,:,:] / 1000. ; # km
        ze2T   = id_mm.variables['e2t'][0,:,:] / 1000. ; # km
                
    (nj,ni) = np.shape(kmaskt)

    kmaskt = np.array(kmaskt, dtype=int)


    zXt = np.zeros((nj,ni))
    zYt = np.zeros((nj,ni))
    zXf = np.zeros((nj,ni))
    zYf = np.zeros((nj,ni))

    # Conversion from Geographic coordinates (lat,lon) to Cartesian in km,
    #  ==> same North-Polar-Stereographic projection as RGPS data...
    zlonT = np.mod( zlonT, 360. )
    zlonF = np.mod( zlonF, 360. )
    zYt[:,:], zXt[:,:] = ConvertGeo2CartesianNPSkm(zlatT, zlonT)
    zYf[:,:], zXf[:,:] = ConvertGeo2CartesianNPSkm(zlatF, zlonF)
    del zlatF, zlonF
    
    # Local resolution in km (for ):
    zResKM = np.zeros((nj,ni))
    zResKM[:,:] = np.sqrt( ze1T*ze1T + ze2T*ze2T )
    del ze1T, ze2T
    #ii = dump_2d_field( 'res_km.nc', zResKM, xlon=zlonT, xlat=zlatT, name='resolution', unit='km' )

    return kmaskt, zlatT, zlonT, zYt, zXt, zYf, zXf, zResKM




def ncSaveCloudBuoys( cf_out, ptime, pIDs, pY, pX, pLat, pLon, mask=[], tunits=tunits_default, fillVal=None ):
    '''
    '''
    print('\n *** [ncSaveCloudBuoys]: About to generate file: '+cf_out+' ...')
    (Nt,) = np.shape(ptime)
    (Nb,) = np.shape(pIDs)
    if np.shape(pY)!=(Nt,Nb) or np.shape(pX)!=(Nt,Nb) or np.shape(pLat)!=(Nt,Nb) or np.shape(pLon)!=(Nt,Nb):
        print('ERROR [ncSaveCloudBuoys]: one of the 2D arrays has a wrong shape!!!')
        exit(0)
    lSaveMask = (np.shape(mask) == (Nt,Nb))
    
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
    v_latb  = f_out.createVariable('latitude' , 'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    v_lonb  = f_out.createVariable('longitude', 'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    v_y     = f_out.createVariable('y_pos' ,    'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    v_x     = f_out.createVariable('x_pos',     'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    if lSaveMask:
        v_mask = f_out.createVariable('mask',   'i1',(cd_time,cd_buoy,),                      zlib=True, complevel=9)

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
        if lSaveMask:
            v_mask[jt,:] = mask[jt,:]
            
    f_out.About  = 'RGPS sea-ice drift data (Kwok, 1998)'
    f_out.Author = 'Generated with `'+path.basename(argv[0])+'` of `mojito` (L. Brodeau, 2022)'
    f_out.close()
    print('      ===> '+cf_out+' saved!')

    return 0


def LoadDist2CoastNC( cNCfile, iverbose=0 ):
    #
    if iverbose>0:
        print('  +++ [util.LoadDist2CoastNC()] Loading "distance to coast" from file:')
        print('       '+cNCfile)
    chck4f(cNCfile)
    with Dataset(cNCfile) as id_in:
        vlon  = id_in.variables['lon'][:]
        vlat  = id_in.variables['lat'][:]
        xdist = id_in.variables['dist'][:,:]
    if iverbose>0: print('       => ok!\n')
        #
    return vlon, vlat, xdist



def LoadNCdataMJT( cfile, krec=0, lmask=False, iverbose=0 ):
    '''
       Open & inspect a `mojito-generated` NetCDF file (like that generated by `ncSaveCloudBuoys`)
       and load the data
    '''
    #
    chck4f(cfile)
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        for cd in ['time','buoy']:
            if not cd in list_dim:
                print(' ERROR [LoadNCdataMJT()]: no dimensions `'+cd+'` found into input file!'); exit(0)

        list_var = list( id_in.variables.keys() )
        for cv in ['id_buoy','latitude','longitude','y_pos','x_pos']:
            if not cv in list_var:
                print(' ERROR [LoadNCdataMJT()]: no variable `'+cv+'` found into input file!'); exit(0)

        Nt = id_in.dimensions['time'].size
        nP = id_in.dimensions['buoy'].size
        if iverbose>0:
            print('\n *** Total number of records in the file = ', Nt)
            print('  *** Total number of buoys in the file = ', nP)
                
        # Time record:
        ctunits = id_in.variables['time'].units
        if not ctunits == tunits_default:
            print(" ERROR [LoadNCdataMJT()]: we expect '"+tunits_default+
                  "' as units for the time record vector, yet we have: "+ctunits)
            exit(0)
        ztime = id_in.variables['time'][krec]

        # Buoys' IDs:
        kBIDs    = np.zeros(nP, dtype=int)
        kBIDs[:] = id_in.variables['id_buoy'][:]
        
        # Coordinates:
        zlat  = id_in.variables['latitude'][krec,:]
        zlon  = id_in.variables['longitude'][krec,:]
        zy    = id_in.variables['y_pos'][krec,:]
        zx    = id_in.variables['x_pos'][krec,:]
        if lmask:
            zmsk = id_in.variables['mask'][krec,:]

    zlon[:] = np.mod(zlon, 360.) ; # Longitudes in the [0:360] frame...

    if lmask:
        return ztime, kBIDs, np.array([zlat,zlon]).T, np.array([zy,zx]).T, zmsk
    else:
        return ztime, kBIDs, np.array([zlat,zlon]).T, np.array([zy,zx]).T


def GetDimNCdataMJT( cfile ):
    '''
    '''
    chck4f(cfile)
    with Dataset(cfile) as id_in:
        list_dim = list( id_in.dimensions.keys() )
        for cd in ['time','buoy']:
            if not cd in list_dim:
                print(' ERROR [GetDimNCdataMJT()]: no dimensions `'+cd+'` found into input file!'); exit(0)
        Nt = id_in.dimensions['time'].size
        nP = id_in.dimensions['buoy'].size
        print(' * [GetDimNCdataMJT]: total number of records in file '+cfile+' =>', Nt)
        print(' * [GetDimNCdataMJT]: max. number of buoys (at start) in file '+cfile+' =>', nP)
    return Nt, nP


