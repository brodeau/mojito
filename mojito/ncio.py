
from sys import argv, exit
from os import path ; #, environ, mkdir
import numpy as np

#from re import split
from netCDF4 import Dataset

from .util import chck4f, epoch2clock, ConvertGeo2CartesianNPSkm



tunits_default = 'seconds since 1970-01-01 00:00:00'

list_required_var_RGPS = [ 'index', 'x', 'y', 'lon', 'lat', 'time', 'idstream', 'streams' ]


def SaveRGPStoNC( ncFile, ptime, pIDs, pykm, pxkm, plat, plon, pqual=[], Nstrm=None, pSid=[] ):

    (nP,) = np.shape(ptime)

    ldoQual   = ( np.shape(pqual)==(nP,) )
    ldoStream = ( Nstrm and np.shape(pSid)==(nP,) )

    with Dataset(ncFile, 'w') as ds:

        points  = ds.createDimension( "points", nP )
        if ldoStream:
            streams = ds.createDimension( "streams", Nstrm )

        time   = ds.createVariable("time", "f8", ("points",), zlib=True, complevel=9)
        x      = ds.createVariable("x", "f8", ("points",), zlib=True, complevel=9)
        y      = ds.createVariable("y", "f8", ("points",), zlib=True, complevel=9)
        lon    = ds.createVariable("lon", "f8", ("points",), zlib=True, complevel=9)
        lat    = ds.createVariable("lat", "f8", ("points",), zlib=True, complevel=9)
        index  = ds.createVariable("index", "long", ("points",), zlib=True, complevel=9)
        if ldoQual:
            qual   = ds.createVariable("q_flag", "byte", ("points",), zlib=True, complevel=9)
        if ldoStream:
            idstream = ds.createVariable("idstream", "byte", ("points",), zlib=True, complevel=9)
            streams = ds.createVariable("streams", "byte", ("streams",) )

        time[:] = ptime[:]
        time.setncattr('units', "seconds since 1970-01-01 00:00:00")
        time.setncattr('standard_name', "time")
        time.setncattr('long_name', "time of virtual buoy")
        time.setncattr('calendar', "standard")

        x[:] = pxkm[:]
        x.setncattr('units', "km")
        x.setncattr('long_name', "X-coordinate in Polar Stereographic projection, lon_0=-45, lat_ts=70")
        x.setncattr('standard_name', "x_coordinate")

        y[:] = pykm[:]
        y.setncattr('units', "km")
        y.setncattr('long_name', "Y-coordinate in Polar Stereographic projection, lon_0=-45, lat_ts=70")
        y.setncattr('standard_name', "y_coordinate")

        lon[:] = plon[:]
        lon.setncattr('units', "degrees_easth")
        lon.setncattr('long_name', "longitude coordinate of virtual buoy")
        lon.setncattr('standard_name', "longitude")

        lat[:] = plat[:]
        lat.setncattr('units', "degrees_north")
        lat.setncattr('long_name', "latitude coordinate of virtual buoy")
        lat.setncattr('standard_name', "latitude")

        index[:] = pIDs[:]
        index.setncattr('long_name', "Id of virtual buoy")

        if ldoQual:
            qual[:] = pqual[:]
            qual.setncattr('long_name', "Quality flag of virtual buoy")

        if ldoStream:
            idstream[:] = pSid[:]
            idstream.setncattr('long_name', "ID of the stream the virtual buoy belongs to")

            streams[:] = np.arange(Nstrm) + 1
            streams.setncattr('long_name', "IDs of existing streams")

        ds.setncattr('title', 'RGPS trajectories')
        ds.setncattr('reference', 'Kwok, Ronald. “The RADARSAT Geophysical Processor System.” (1998).')

    print(' *** File '+ncFile+' generated!\n')

    return 0




def LoadDataRGPS( cfile, iverbose=0 ):
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
        for cv in list_required_var_RGPS:
            if not cv in list_var:
                print(' ERROR [LoadDataRGPS()]: no variable `'+cv+'` found into input file!'); exit(0)

        nP = id_in.dimensions['points'].size
        if iverbose>0: print('\n *** Total number of points in the file = ', nP)
        nS = id_in.dimensions['streams'].size
        if iverbose>0: print('\n *** Total number of streams in the file = ', nS)

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

        # Buoy stream:
        kStrm    = np.zeros(nP, dtype='i1')
        kStrm[:] = id_in.variables['idstream'][:]

    zlon[:] = np.mod(zlon, 360.) ; # Longitudes in the [0:360] frame...

    return nP, nS, ztime, zy, zx, zlat, zlon, kBIDs, kStrm



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


def GetModelUVGrid( fNCmeshmask ):

    chck4f( fNCmeshmask)

    # Reading mesh metrics into mesh-mask file:
    with Dataset(fNCmeshmask) as id_mm:
        zlonV = id_mm.variables['glamv'][0,:,:]
        zlatV = id_mm.variables['gphiv'][0,:,:]
        zlonU = id_mm.variables['glamu'][0,:,:]
        zlatU = id_mm.variables['gphiu'][0,:,:]

    (nj,ni) = np.shape(zlonV)

    zXv = np.zeros((nj,ni))
    zYv = np.zeros((nj,ni))
    zXu = np.zeros((nj,ni))
    zYu = np.zeros((nj,ni))

    # Conversion from Geographic coordinates (lat,lon) to Cartesian in km,
    #  ==> same North-Polar-Stereographic projection as RGPS data...
    zlonV = np.mod( zlonV, 360. )
    zlonU = np.mod( zlonU, 360. )
    zYv[:,:], zXv[:,:] = ConvertGeo2CartesianNPSkm(zlatT, zlonV)
    zYf[:,:], zXu[:,:] = ConvertGeo2CartesianNPSkm(zlatF, zlonU)
    del zlatF, zlonU

    return zYv, zXv, zYu, zXu




def ncSaveCloudBuoys( cf_out, ptime, pIDs, pY, pX, pLat, pLon, mask=[], xtime=[],
                      tunits=tunits_default, fillVal=None, corigin=None ):
    '''
    '''
    print('\n *** [ncSaveCloudBuoys]: About to generate file: '+cf_out+' ...')
    (Nt,) = np.shape(ptime)
    (Nb,) = np.shape(pIDs)
    if np.shape(pY)!=(Nt,Nb) or np.shape(pX)!=(Nt,Nb) or np.shape(pLat)!=(Nt,Nb) or np.shape(pLon)!=(Nt,Nb):
        print('ERROR [ncSaveCloudBuoys]: one of the 2D arrays has a wrong shape!!!')
        exit(0)
    lSaveMask = (np.shape(mask)  == (Nt,Nb))
    lSaveTime = (np.shape(xtime) == (Nt,Nb)) ; # time for each buoy!
    #
    f_out = Dataset(cf_out, 'w', format='NETCDF4')
    #
    # Dimensions:
    cd_time = 'time'
    cd_buoy = 'buoy'
    f_out.createDimension(cd_time, None)
    f_out.createDimension(cd_buoy, Nb  )
    #
    # Variables:
    v_time = f_out.createVariable(cd_time,     'i4',(cd_time,))
    v_buoy = f_out.createVariable(cd_buoy,     'i4',(cd_buoy,))
    v_bid  = f_out.createVariable('id_buoy',   'i4',(cd_buoy,))
    x_lat  = f_out.createVariable('latitude' , 'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    x_lon  = f_out.createVariable('longitude', 'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    x_ykm  = f_out.createVariable('y_pos' ,    'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    x_xkm  = f_out.createVariable('x_pos',     'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    #
    v_time.units = tunits
    v_bid.units  = 'ID of buoy'
    x_lat.units = 'degrees north'
    x_lon.units = 'degrees south'
    x_ykm.units    = 'km'
    x_xkm.units    = 'km'
    #
    if lSaveMask:
        v_mask = f_out.createVariable('mask',  'i1',(cd_time,cd_buoy,),                      zlib=True, complevel=9)
    if lSaveTime:
        x_tim  = f_out.createVariable('time_pos', 'i4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
        x_tim.units = tunits
    #
    v_buoy[:] = np.arange(Nb,dtype='i4')
    v_bid[:]  = pIDs[:]
    #
    for jt in range(Nt):
        v_time[jt]  = ptime[jt]
        x_lat[jt,:] =  pLat[jt,:]
        x_lon[jt,:] =  pLon[jt,:]
        x_ykm[jt,:] =    pY[jt,:]
        x_xkm[jt,:] =    pX[jt,:]
        if lSaveMask:
            v_mask[jt,:] = mask[jt,:]
        if lSaveTime:
            x_tim[jt,:] = xtime[jt,:]
    #
    if corigin:
         f_out.Origin = corigin
    f_out.About  = 'Lagrangian sea-ice drift'
    f_out.Author = 'Generated with `'+path.basename(argv[0])+'` of `mojito` (L. Brodeau, 2023)'
    f_out.close()
    print('      ===> '+cf_out+' saved!')
    #
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



def LoadNCtimeMJT( cfile, iverbose=0 ):
    '''
       Open & inspect a `mojito-generated` NetCDF file (like that generated by `ncSaveCloudBuoys`)
       and load the time vector only !
    '''
    #
    chck4f(cfile)
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        if not 'time' in list_dim:
            print(' ERROR [LoadNCtimeMJT()]: no dimensions `'+time+'` found into input file!'); exit(0)
        #
        list_var = list( id_in.variables.keys() )
        if not 'time' in list_var:
            print(' ERROR [LoadNCtimeMJT()]: no variable `'+time+'` found into input file!'); exit(0)

        Nt = id_in.dimensions['time'].size

        # Time record:
        ctunits = id_in.variables['time'].units
        if not ctunits == tunits_default:
            print(" ERROR [LoadNCtimeMJT()]: we expect '"+tunits_default+
                  "' as units for the time record vector, yet we have: "+ctunits)
            exit(0)
        ztime = id_in.variables['time'][:]

        print(' * [LoadNCtimeMJT()]: read the "time" of the '+str(Nt)+' records in file '+path.basename(cfile))

    return Nt, ztime



def LoadNCdataMJT( cfile, krec=0, lmask=False, lGetTimePos=False, iverbose=0 ):
    '''
       Open & inspect a `mojito-generated` NetCDF file (like that generated by `ncSaveCloudBuoys`)
       and load the data
    '''
    #
    list_var_needed = ['id_buoy','latitude','longitude','y_pos','x_pos']
    if lGetTimePos:
        list_var_needed.append('time_pos')

    chck4f(cfile)
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        for cd in ['time','buoy']:
            if not cd in list_dim:
                print(' ERROR [LoadNCdataMJT()]: no dimensions `'+cd+'` found into input file!'); exit(0)

        list_var = list( id_in.variables.keys() )
        for cv in list_var_needed:
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
        if lGetTimePos:
            ztpos = id_in.variables['time_pos'][krec,:]

    zlon[:] = np.mod(zlon, 360.) ; # Longitudes in the [0:360] frame...

    if lmask:
        if lGetTimePos:
            return ztime, kBIDs, np.array([zlat,zlon]).T, np.array([zy,zx]).T, zmsk, ztpos
        else:
            return ztime, kBIDs, np.array([zlat,zlon]).T, np.array([zy,zx]).T, zmsk
    else:
        if lGetTimePos:
            return ztime, kBIDs, np.array([zlat,zlon]).T, np.array([zy,zx]).T, ztpos
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
        corgn = id_in.Origin
        print('   * [GetDimNCdataMJT]: total number of records in file '+cfile+' =>', Nt)
        print('   * [GetDimNCdataMJT]: max. number of buoys (at start) in file '+cfile+' =>', nP)
        print('   * [GetDimNCdataMJT]: origin of data: "'+corgn+'"')
        list_var = list( id_in.variables.keys() )
        ltimePos = ('time_pos' in list_var)
    return Nt, nP, corgn, ltimePos



def SeedFileTimeInfo( fSeedNc, iverbose=0 ):
    from re import split
    #
    cSeed = str.replace( path.basename(fSeedNc), 'SELECTION_', '' )
    cSeed = str.replace( cSeed, '.nc', '' )
    cBtch = split('_',path.basename(fSeedNc))[2]

    chck4f(fSeedNc)
    print('\n *** Will read initial seeding positions in first record of file:\n      => '+fSeedNc+' !')
    ntr, zt = LoadNCtimeMJT( fSeedNc, iverbose=iverbose )
    idate0 = zt[0]     ; cdate0 = epoch2clock(idate0)
    idateN = zt[ntr-1] ; cdateN = epoch2clock(idateN)
    print('    => earliest and latest time position in the file: '+cdate0+' - '+cdateN)

    zdDate =  int( round( (idateN - idate0)/3600., 0 ) * 3600. )
    print('    => rounded time span =>',zdDate/3600.,'hours')
    idate0 =  int( round( idate0/3600., 0 ) * 3600. ) ; cdate0 = epoch2clock(idate0)
    idateN = int( idate0 + zdDate )  ; cdateN = epoch2clock(idateN)
    print('    ==> will actually use rounded to the hour! => '+cdate0+' - '+cdateN)

    return idate0, idateN, cSeed, cBtch


def ModelFileTimeInfo( fModelNc, iverbose=0 ):
    #
    from re import split
    #
    with Dataset(fModelNc) as ds_mod:
        Nt = ds_mod.dimensions['time_counter'].size
        if ds_mod.variables['time_counter'].units != tunits_default:
            print('ERROR: wrong units for time calendar in file:',fModelNc)
            exit(0)
        ztime = np.array( ds_mod.variables['time_counter'][:] , dtype='i4' )
    #
    print('\n\n *** '+str(Nt)+' records in input MODEL file!')
    #
    idate0, idateN = np.min(ztime), np.max(ztime)
    print('    * [ModelFileTimeInfo] First and last dates in MODEL input file:',epoch2clock(idate0), epoch2clock(idateN))
    #
    # Infer name of NEMO CONFIG and experiment from SI3 file:
    vn = split('_',path.basename(fModelNc))
    nconf, nexpr = vn[0], split('-',vn[1])[1]
    print('    * [ModelFileTimeInfo] NEMO config and experiment =', nconf, nexpr)
    #
    return Nt, ztime, idate0, idateN, nconf, nexpr

