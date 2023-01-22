#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import exit
#from os import path, mkdir
import numpy as np



def LoadDist2CoastNC( cNCfile, ivrbs=0 ):
    from climporn import chck4f
    from netCDF4  import Dataset
    #
    if ivrbs>0:
        print('  +++ [util.LoadDist2CoastNC()] Loading "distance to coast" from file:')
        print('       '+cNCfile)
    chck4f(cNCfile)
    with Dataset(cNCfile) as id_in:
        vlon  = id_in.variables['lon'][:]
        vlat  = id_in.variables['lat'][:]
        xdist = id_in.variables['dist'][:,:]
    if ivrbs>0: print('       => ok!\n')
        #
    return vlon, vlat, xdist


def Dist2Coast( lon0, lat0, plon, plat, pdist2coat ):
    '''
       Returns the distance to the nearest coast of a given point (lon,lat)
        INPUT:
          * lon0, lat0: coordinates (scalars) [degrees East], [degrees West]
          * plon:       1D array of longitudes assosiated to `pdist2coat`
          * plat:       1D array of latitudes  assosiated to `pdist2coat`
          * pdist2coat: 2D array containing rasterized distance to coast (originally read in NC file) [km]
    '''
    from climporn import degE_to_degWE
    #
    rx = np.mod( lon0, 360. ) ; # [0:360] frame
    vx = np.mod( plon, 360. ) ; # [0:360] frame
    # Are we dangerously close to the [0--360] cut?
    #  => then will work in the [-180:180] frame:
    if rx > 355.:
        rx = degE_to_degWE( rx )
        vx = degE_to_degWE( vx )
        #print(' lon0, lat0 =', rx, lat0)
    ip = np.argmin( np.abs(  vx[:] - rx  ) )
    jp = np.argmin( np.abs(plat[:] - lat0) )
    #print(' ip, jp =', ip, jp)
    #print(' Nearest lon, lat =', plon[ip], plat[jp])
    del vx, rx
    return max( pdist2coat[jp,ip] , 0. )


def TimeBins4Scanning( pdt1, pdt2, pdt,  iverbose=0 ):
    '''
       Need a time axis with bins to look for buoys within...

       Input:
                * pdt1, pdt2 : start & end time ([s] UNIX epoch time)
                * pdt        : width of bin ([s] UNIX epoch time)
       Returns:
                * nB :  number of bins
                * vTB : array(nB,3) axe==0 => time at center of bin      ([s] UNIX epoch time)
                *                        axe==1 => time at lower bound of bin ([s] UNIX epoch time)
                *                        axe==2 => time at upper bound of bin ([s] UNIX epoch time)
                * cTc : array(nB) human-readable date corresponding to time at center of bin (character string)
    '''
    from climporn import epoch2clock
    #
    zhdt = pdt/2.
    #
    nB = int(round((pdt2 - pdt1) / pdt))
    print('\n *** New fixed time axis to use to scan data:\n    ===> nB = '+str(nB)+' time bins!')
    #
    vTB = np.zeros((nB,3), dtype=int  ) ; # `*,0` => precise time | `*,1` => bound below | `*,2` => bound above
    cTc = np.zeros( nB   , dtype='U19')
    #
    vTB[0,0] =  pdt1 + zhdt    ; # time at center of time bin
    cTc[0]   = epoch2clock(vTB[0,0])
    for jt in range(1,nB):
        tt = vTB[jt-1,0] + pdt
        vTB[jt,0] = tt                ; # time at center of time bin
        cTc[jt]   = epoch2clock(tt)
    # Time bins bounds:
    vTB[:,1] = vTB[:,0] - zhdt
    vTB[:,2] = vTB[:,0] + zhdt
    #
    if iverbose>0:
        for jt in range(nB):
            print(" --- jt="+'%3.3i'%(jt)+": * center of bin => ",vTB[jt,0]," => ",cTc[jt])
            print("             * bin bounds    => "+epoch2clock(vTB[jt,1])+" - "+epoch2clock(vTB[jt,2])+"\n")
    #
    return nB, vTB, cTc


def InspectLoadData( cfile, plistVar ):
    '''
       Open & inspect the RGPS NetCDF input file and load the raw data
    '''
    from netCDF4  import Dataset
    from climporn import chck4f
    #
    ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!
    #
    chck4f(cfile)
    #
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        if not 'points' in list_dim:
            print(' ERROR: no dimensions `points` found into input file!'); exit(0)

        list_var = list( id_in.variables.keys() ) ; print(' ==> variables: ', list_var, '\n')
        for cv in plistVar:
            if not cv in list_var:
                print(' ERROR: no variable `'+cv+'` found into input file!'); exit(0)

        nP = id_in.dimensions['points'].size
        print('\n *** Total number of points in the file = ', nP)

        # Time records:
        ctunits = id_in.variables['time'].units
        if not ctunits == ctunits_expected:
            print(" ERROR: we expect '"+ctunits_expected+"' as units for the time record vector, yet we have: "+ctunits)
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


def KeepDataInterest( pdt1, pdt2, ptime0, pBIDs0, px0, py0, plon0, plat0,  rmskVal=-99999. ):
    '''
       Only keep data of interest in 1D arrays (except time array), based on date1 and date2.
       Excluded points are masked...
       Again, time array `ptime0` won't be masked (unnecessary and dangerous)
       Input:
                * pdt1, pdt2 : start & end time ([s] UNIX epoch time)

       Returns:
                * nB   : number of buoys that exist during the specified time range
                * zIDs : array(nB), IDs of the buoys that exist during the specified time range
    '''
    nP = len(pBIDs0)

    # Will mask all point that are before and beyond our period of interest:
    zmsk = np.zeros(nP, dtype=int) + 1
    zmsk[np.where(ptime0 < pdt1)] = 0
    zmsk[np.where(ptime0 > pdt2)] = 0
    #
    (idx_masked,) = np.where( zmsk == 0 )
    #
    if nP-len(idx_masked) != np.sum(zmsk):
        print('ERROR: fuck up #1!')
        exit(0)
    print('\n *** Total number of points remaining after time-range-exclusion = ',nP-len(idx_masked), '=', np.sum(zmsk))
    #
    pBIDs0[idx_masked] = int(rmskVal) ; pBIDs0 =  np.ma.masked_where( zmsk==0, pBIDs0 )
    px0[idx_masked]    =     rmskVal  ; px0    =  np.ma.masked_where( zmsk==0, px0    )
    py0[idx_masked]    =     rmskVal  ; py0    =  np.ma.masked_where( zmsk==0, py0    )
    plon0[idx_masked]  =     rmskVal  ; plon0  =  np.ma.masked_where( zmsk==0, plon0  )
    plat0[idx_masked]  =     rmskVal  ; plat0  =  np.ma.masked_where( zmsk==0, plat0  )
    #
    # Remaining buoys (IDs)
    (idx,) = np.where(pBIDs0.data > 0)
    zIDs = np.sort( np.unique( pBIDs0[idx] ) ) ; # if not `[idx]` then `rmskVal` is counted once!!!
    nB   = len(zIDs)
    print("\n *** We found "+str(nB)+" different buoys alive during specified period of time!")
    #
    return nB, zIDs


def ValidCnsctvRecordsBuoy( kid, ptime0, pBIDs0, dt_expected, max_dev_from_dt_expected ):
    '''
         * Analysis of the time records for this particular buoy...
            => Must cut off the series of the buoys as soon as its dt is too far
               from the nominal time step:
           => based on the initial time for this particular buoy:
              - construct the ideal expected `ztime` (`ztime_ideal`) based on nominal dt
              - dezing tout ce qui s'eloigne trop de ce ztime_ideal !
       Input:
                * kid    : ID of treated buoy
                * ptime0 : time array as in raw data
                * pBIDs0 : IDs array as in raw data
                * dt_expected: expected time gap between 2 consecutive records of a buoy [s]
                * max_dev_from_dt_expected: overshoot tolerance for `dt_expected` [s]
       Returns:
               * nbROK   : number of valid consecutive records for this buoy
               * idx0_id : array of location indices (in the raw data arrays) for these valid records of this buoy
               * ztime   : array of dates associated with all these records [s]
    '''
    (idx_id,) = np.where( pBIDs0 == kid)
    #
    ztime  = ptime0[idx_id] ; # all time records for this particular buoy
    nbR1b = len(ztime)      ; # n. of time records for this particulat buoy
    #
    nbROK = nbR1b
    ztime_ideal = np.array( [ ztime[0]+float(i)*float(dt_expected) for i in range(nbR1b) ], dtype=float )
    vtdev = np.abs(ztime - ztime_ideal)
    if np.any(vtdev > max_dev_from_dt_expected):
        (indFU,) = np.where(vtdev > max_dev_from_dt_expected)
        nbROK = np.min(indFU) ; # yes! no -1 !!!
    #
    if nbROK < nbR1b:
        # Update with only keeping acceptable time records (#fixme: for now only those until first fuckup)
        idx_id = idx_id[0:nbROK]
        ztime   =   ztime[0:nbROK]
    #
    del nbR1b, ztime_ideal, vtdev
    #
    if len(idx_id)!=nbROK or len(ztime)!=nbROK:
        print('ERROR: ValidCnsctvRecordsBuoy => `len(idx_id)!=nbROK or len(ztime)!=nbROK`')
        exit(0)
    #
    return nbROK, idx_id, ztime



def OrderCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''

    # sort the points based on their x-coordinates
    isortX  = np.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]

    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]

    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = np.argsort(leftMost[:,1])
    leftMost   = leftMost[isortL,:]
    (tl, bl)   = leftMost
    [i1l, i2l] = isortML[isortL]

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = np.argsort(rghtMost[:,1])
    rghtMost   = rghtMost[isortR,:]
    (tr, br)   = rghtMost
    [i1r, i2r] = isortMR[isortR]

    #idx = np.concatenate([idxL,idxR])
    #isort = np.array([isortX[i] for i in np.concatenate([isortL,isortR])])

    del isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR

    # return the coordinates in top-left, top-right,
    # bottom-right, and bottom-left order
    return np.array([tl, tr, br, bl], dtype="float32") , np.array([i1l, i1r, i2r, i2l])


def OrderCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''
    zz = xcoor.copy()
    iz = np.array([0,1,2,3])
    zt, it = OrderCW(xcoor)
    #
    zz[1:4,:] = zt[:0:-1,:]
    iz[1:4]   = it[:0:-1]
    del zt, it
    return zz, iz



def SortIndicesCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter-clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''
    # sort the points based on their x-coordinates
    isortX  = np.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]

    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]

    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = np.argsort(leftMost[:,1])
    [i1l, i2l] = isortML[isortL]

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = np.argsort(rghtMost[:,1])
    [i1r, i2r] = isortMR[isortR]

    #idx = np.concatenate([idxL,idxR])
    #isort = np.array([isortX[i] for i in np.concatenate([isortL,isortR])])

    del xSorted, isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR

    return np.array([i1l, i1r, i2r, i2l])




def mergeNPZ( list_npz_files, cf_out='merged_file.npz', iverbose=0 ):
    '''
       Merge several npz files of the same date into a single one!
    '''
    from climporn import epoch2clock

    nbf = len(list_npz_files)

    # First round to check the number of points in each file
    vit, vNbP = [], []
    for cf in list_npz_files:
        with np.load(cf) as data:
            it = int( data['itime'] ) ; # because otherwize it's an array of size 1 ?
            vit.append(it)
            npf = len(data['vids'])
            vNbP.append(npf)
        if iverbose>0: print('  '+cf+' => time ='+epoch2clock(it)+' | '+str(npf)+' points!')
    #
    # For merged file:
    it = int( round(np.mean(vit),0) )
    nP = np.sum(vNbP)
    if iverbose>0: print(' vNbP =', vNbP,'=>',nP,'points in total!, time =',epoch2clock(it))

    vx, vy = np.zeros(nP), np.zeros(nP)
    vlon, vlat = np.zeros(nP), np.zeros(nP)
    vids   = np.zeros(nP, dtype=int)

    jf=0
    i1, i2 = 0, 0
    for cf in list_npz_files:
        npf = vNbP[jf]
        i2=i2+npf
        with np.load(cf) as data:
            vx[i1:i2], vy[i1:i2] = data['vx'], data['vy']
            vlon[i1:i2], vlat[i1:i2] = data['vlon'], data['vlat']
            vids[i1:i2]          = data['vids']
        #
        i1=i1+npf
        jf=jf+1

    if len(np.unique(vids)) < nP:
        print(' ERROR: [util.mergeNPZ] => some IDs are identical :('); exit(0)

    # Time to save in the new npz file:
    if iverbose>0: print('  [util.mergeNPZ] ==> saving merged files into '+cf_out+' !')
    np.savez_compressed( cf_out, itime=it, vx=vx, vy=vy, vlon=vlon, vlat=vlat, vids=vids )

    return 0


def idx_suppress_xy_copies( X, rmask_val=-999. ):
    '''
         Input: `X` is a 2D array of float containing x,y coordinates for Np points => shape=(Np,2)

             ==> will return the 1D integer vector containing row (points) indices (along Np)
                 corresponding to the location of points that already exist once, so the can be
                 suppresses or masked later on... (only a single occurence of a point coord. is
                 kept). While ignoring masked values flagged with `rmask_val`
    '''
    _,nc = np.shape(X)
    if (nc !=2 ):
        print('ERROR [idx_suppress_xy_copies]: second dimmension of X must be 2! (coordinates)')

    # Row indices that exclude masked points:
    (Ix,_) = np.where( X != [rmask_val,rmask_val] )
    Ix = Ix[::2]

    # Row indices that select masked points:
    (Iz,_) = np.where( X == [rmask_val,rmask_val] )
    Iz = Iz[::2]

    _,Iu = np.unique( X, axis=0, return_index=True )

    # keep the values of `Iu` that are not in `Iz`:
    Ic = np.setdiff1d( Iu, Iz ) ;

    # keep the values of `Ix` that are not in `Ic`:
    return np.setdiff1d( Ix, Ic )



def SubSampCloud( rd_km, pCoor, pIDs, ptime, pNames=[] ):
    '''
    '''
    from gudhi import subsampling as sbspl
    #
    cerr = 'ERROR [util.SubSampCloud()]: '
    l_do_names = ( len(pNames) > 0 )

    # Sanity check of input:
    if rd_km <= 0. or rd_km > 2000:
        print(cerr+'silly value for `rd_km`:',rd_km)
        exit(0)
    (Nb0,n2) = np.shape(pCoor)
    if (n2 != 2 ):
        print(cerr+'second dimmension of `pCoor` must be 2 !')
        exit(0)
    if len(pIDs) != Nb0:
        print(cerr+'len(pIDs) != Nb0 !')
        exit(0)
    if len(ptime) != Nb0:
        print(cerr+'len(ptime) != Nb0 !')
        exit(0)
    if l_do_names:
        if len(pNames) != Nb0:
            print(cerr+'len(pNames) != Nb0 !')
            exit(0)

    zCoor = np.array( sbspl.sparsify_point_set( pCoor, min_squared_dist=rd_km*rd_km ) )
    (Nb,_) = np.shape(zCoor)

    # Retrieve corresponding indices for selected points:
    ileft = np.zeros(Nb, dtype=int)
    for i in range(Nb):
        (idx,_) = np.where( pCoor[:,:]==zCoor[i,:] )
        ileft[i] = idx[0]

    if l_do_names:
        return Nb, zCoor, pIDs[ileft], ptime[ileft], pNames[ileft]
    else:
        return Nb, zCoor, pIDs[ileft], ptime[ileft]




def StdDev( pmean, pX ):
    zz = pX[:] - pmean
    return np.sqrt( np.mean( zz*zz ) )




