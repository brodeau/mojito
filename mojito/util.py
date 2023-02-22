#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import exit
#from os import path ; #, mkdir
import numpy as np


def chck4f( cf ):
    from os.path import exists
    if not exists(cf):
        print(' ERROR [chck4f()]: file '+cf+' does not exist!')
        exit(0)

def degE_to_degWE( X ):
    '''
    # From longitude in 0 -- 360 frame to -180 -- +180 frame...
    '''
    if np.shape( X ) == ():
        # X is a scalar
        from math import copysign
        return     copysign(1., 180.-X)*        min(X,     abs(X-360.))
    else:
        # X is an array
        return np.copysign(1., 180.-X)*np.minimum(X, np.abs(X-360.))



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


def KeepDataInterest( pdt1, pdt2, ptime0, pBIDs0, px0, py0, plon0, plat0,  rmskVal=-99999. ):
    '''
       Only keep data of interest in 1D arrays (except time array), based on date1 and date2.
       Excluded points are masked...
       Again, time array `ptime0` won't be masked (unnecessary and dangerous)
       Input:
                * pdt1, pdt2 : start & end time ([s] UNIX epoch time)

       Returns:
                * nB   : number of different buoys that exist for at least 1 record during specified date range
                * zIDs : list (unique) of IDs for these buoys (array[nB] of int)
    '''
    nP = len(pBIDs0)
    #
    # Will mask all point that are before and beyond our period of interest:
    zmsk = np.zeros(nP, dtype=int) + 1
    zmsk[np.where(ptime0 < pdt1)] = 0
    zmsk[np.where(ptime0 > pdt2)] = 0
    #
    (idx_masked,) = np.where( zmsk == 0 )
    #
    if nP-len(idx_masked) != np.sum(zmsk):
        print('ERROR: [KDI] fuck up #1!')
        exit(0)
    print('\n *** [KDI] Total number of points remaining after time-range-exclusion = ',nP-len(idx_masked), '=', np.sum(zmsk))
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
    print("\n *** [KDI] We found "+str(nB)+" different buoys alive during specified period of time!")
    #
    return nB, zIDs


def ValidCnsctvRecordsBuoy( time_min, kidx, ptime0, pBIDs0, pidx_ignore, dt_expected, max_dev_from_dt_expected ):
    '''
    LOLO: add test to see if time vector we return is increasing!!! and nothing stupid / boundaries!!!

         * Analysis of the time records for this particular buoy...
            => Must cut off the series of the buoys as soon as its dt is too far
               from the nominal time step:
           => based on the initial time for this particular buoy:
              - construct the ideal expected `ztime` (`ztime_ideal`) based on nominal dt
              - dezing tout ce qui s'eloigne trop de ce ztime_ideal !
       Input:
                * time_min: all buoys considered should exist at and after this time, not before!!! epoch time in [s]
                * kidx   : the indice of the point we are dealing with (in the frame of original raw data)
                * ptime0 : time array as in raw data
                * pBIDs0 : IDs array as in raw data
                * dt_expected: expected time gap between 2 consecutive records of a buoy [s]
                * max_dev_from_dt_expected: overshoot tolerance for `dt_expected` [s]
       Returns:
               * nbROK   : number of valid consecutive records for this buoy
               * idx0_id : array of location indices (in the raw data arrays) for these valid records of this buoy
               * ztime   : array of dates associated with all these records [s]
    '''
    from climporn import epoch2clock
    #
    (idx_buoy,   ) = np.where( pBIDs0 == pBIDs0[kidx] )
    (idx_exclude,) = np.where( ptime0 < time_min )
    idx_keep0      = np.setdiff1d( idx_buoy,  idx_exclude ) ; # keep the values of `idx_buoy` that are not in `idx_exclude`
    idx_keep       = np.setdiff1d( idx_keep0, pidx_ignore ) ; # keep the values of `idx_keep0` that are not in `pidx_ignore`
    del idx_exclude, idx_buoy, idx_keep0
    #
    ztime  = ptime0[idx_keep] ; # all time records for this particular buoy
    if np.any( ztime<=time_min ):
        print('ERROR: ValidCnsctvRecordsBuoy => `np.any( ztime<=time_min )` !!!!'); exit(0)
    #
    nbR1b = len(idx_keep)      ; # n. of time records for this particulat buoy
    #
    nbROK = nbR1b
    #
    # We cut the series at the first occurence of a dt (time between 2 consec. pos.) > `max_dev_from_dt_expected`
    ztime_ideal = np.array( [ ztime[0]+float(i)*float(dt_expected) for i in range(nbR1b) ], dtype=float )
    vtdev = np.abs(ztime - ztime_ideal)
    if np.any(vtdev > max_dev_from_dt_expected):
        #print('LOLO: vtdev[:]/3600. = ',vtdev[:]/3600.)
        (indFU,) = np.where(vtdev > max_dev_from_dt_expected)
        nbROK = np.min(indFU) ; # yes! no -1 !!!
        #print('LOLO => nbROK =',nbROK)
    #
    if nbROK < nbR1b:
        # Update with only keeping acceptable time records (#fixme: for now only those until first fuckup)
        idx_keep = idx_keep[0:nbROK]
        ztime   =   ztime[0:nbROK]
    #
    return nbROK, idx_keep, ztime



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




def mergeNPZ( list_npz_files, t_ref, cf_out='merged_file.npz', iverbose=0 ):
    '''
       Merge several npz files of the "considered identical" date into a single one!
        * t_ref: reference time for these files (epoch UNIX) [s]
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
    itime_mean = int( round(np.mean(vit),0) )
    nP = np.sum(vNbP)
    if iverbose>0: print(' vNbP =', vNbP,'=>',nP,'points in total!, time =',epoch2clock(itime_mean))

    vtime  = np.zeros(nP)
    vx, vy = np.zeros(nP), np.zeros(nP)
    vlon, vlat = np.zeros(nP), np.zeros(nP)
    vids   = np.zeros(nP, dtype=int)

    nPmax = np.max(vNbP)
    xids   = np.zeros((nbf,nPmax),dtype=int) - 999

    jf=0
    i1, i2 = 0, 0
    for cf in list_npz_files:
        npf = vNbP[jf]
        i2=i2+npf
        with np.load(cf) as data:
            vtime[i1:i2]         = data['vtime']
            vx[i1:i2], vy[i1:i2] = data['vx'], data['vy']
            vlon[i1:i2], vlat[i1:i2] = data['vlon'], data['vlat']
            vids[i1:i2]          = data['vids']
        #
        xids[jf,:npf] = vids[i1:i2]
        #
        i1=i1+npf
        jf=jf+1

    # There might be IDs present more than once (that should not happen if initial scrips was good?)
    vids_u, idxvu = np.unique(vids, return_index=True )

    if len(vids_u) < nP:
        print(' WARNING: [util.mergeNPZ] => some IDs are identical :(',nP-len(vids_u),'of them...')

        vids_ref = vids.copy()
        vidx     = np.arange(nP)

        nP, zidx = SuppressMulitOccurences( vids, vtime, vids_ref, vidx, t_ref )

        vids  = vids_ref[zidx]
        vtime = vtime[zidx]
        vx    = vx[zidx]
        vy    = vy[zidx]
        vlon  = vlon[zidx]
        vlat  = vlat[zidx]

        del vids_ref, vidx, zidx

    if epoch2clock(itime_mean) != epoch2clock(t_ref):
        print(' WARNING: [util.mergeNPZ] => `epoch2clock(itime_mean) != epoch2clock(t_ref)` !!!')

    # Time to save in the new npz file:
    if iverbose>0: print('  [util.mergeNPZ] ==> saving merged files into '+cf_out+' !')
    np.savez_compressed( cf_out, itime=itime_mean, date=epoch2clock(itime_mean), Npoints=nP, vids=vids,
                         vtime=vtime, vx=vx, vy=vy, vlon=vlon, vlat=vlat,  )
    return 0


def IndXYClones( pXY, rmask_val=-999. ):
    '''
         Input: `pXY` is a 2D array of float containing x,y coordinates for Np points => shape=(Np,2)

             ==> will return the 1D integer vector containing row (points) indices (along Np)
                 corresponding to the location of points that already exist once, so the can be
                 suppresses or masked later on... (only a single occurence of a point coord. is
                 kept). While ignoring masked values flagged with `rmask_val`
    '''
    (_,nc) = np.shape(pXY)
    if (nc !=2 ):
        print('ERROR [IndXYClones]: second dimmension of pXY must be 2! (coordinates)')
        exit(0)

    # Row indices that exclude masked points:
    (IXvalid,_) = np.where( pXY != [rmask_val,rmask_val] )
    IXvalid = IXvalid[::2]

    # Row indices that select masked points:
    (IXmaskd,_) = np.where( pXY == [rmask_val,rmask_val] )
    IXmaskd = IXmaskd[::2]

    _,Iunq = np.unique( pXY, axis=0, return_index=True )

    # keep the values of `Iunq` that are not in `IXmaskd`:
    IXokunq = np.setdiff1d( Iunq, IXmaskd ) ;

    # keep the values of `IXvalid` that are not in `IXokunq`:
    # => that's the indices of points that could be removed to shrink `pXY` !
    return np.setdiff1d( IXvalid, IXokunq )


def GetRidOfXYClones( pXY, rmask_val=-999. ):
    '''
         Input: `pXY` is a 2D array of float containing x,y coordinates for Np points => shape=(Np,2)

             ==> will return updated version of coordinates array `pXY`, rid of points already existing (clones)
    '''
    (_,nc) = np.shape(pXY)
    if (nc !=2 ):
        print('ERROR [GetRidOfXYClones]: second dimmension of pXY must be 2! (coordinates)')
        exit(0)

    # Row indices that select masked points:
    (IXmaskd,_) = np.where( pXY == [rmask_val,rmask_val] )
    IXmaskd = IXmaskd[::2]

    _,Iunq = np.unique( pXY, axis=0, return_index=True )

    # keep the values of `Iunq` that are not in `IXmaskd`:
    IXokunq = np.setdiff1d( Iunq, IXmaskd ) ;

    # keep the values of `IXvalid` that are not in `IXokunq`:
    # => that's the indices of points that could be removed to shrink `pXY` !
    return pXY()



def SubSampCloud( rd_km, pCoor ):
    '''
       * pCoor: [X,Y] !!!
    '''
    from gudhi import subsampling as sbspl
    #
    cerr = 'ERROR [util.SubSampCloud()]: '
    # Sanity check of input:
    if rd_km <= 0. or rd_km > 2000:
        print(cerr+'silly value for `rd_km`:',rd_km)
        exit(0)
    (Nb0,n2) = np.shape(pCoor)
    if (n2 != 2 ):
        print(cerr+'second dimmension of `pCoor` must be 2 !')
        exit(0)
    
    zCoor = np.array( sbspl.sparsify_point_set( pCoor, min_squared_dist=rd_km*rd_km ) )
    (Nb,_) = np.shape(zCoor)

    # Retrieve corresponding indices for selected points:
    idxleft = np.zeros(Nb, dtype=int)
    for i in range(Nb):
        (idx,_) = np.where( pCoor[:,:]==zCoor[i,:] )
        idxleft[i] = idx[0]

    return Nb, zCoor, idxleft
    #if l_do_cgeo:
    #    zgeo = np.array( [ pLonLat[ileft,0], pLonLat[ileft,1] ] ).T
    #    if l_do_name:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft], zgeo, pNames[ileft]
    #    else:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft], zgeo
    #else:
    #    if l_do_name:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft], pNames[ileft]
    #    else:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft]




def StdDev( pmean, pX ):
    zz = pX[:] - pmean
    return np.sqrt( np.mean( zz*zz ) )



def SuppressMulitOccurences( pIDs, ptime, pIDsRef0, pidx0, rtime, iverbose=0 ):
    '''
         For many possible reasons the same buoy ID can exist more than once in `pIDs`,
         => we need to keep only one occurence of the location index of these points,
            based on the date (closest to center of bin `rTc`)
    INPUT:
            * pIDs     : 1D array of integers containing buoy IDs with muli-occurence of certain IDs
            * ptime    : 1D array of real containing epoch time date [s] associated to each buoy
            * pIDsRef0 : 1D array of integers containing buoy IDs the "0" reference
            * pidx0     : 1D array of integers containing indices that do this: pIDs == pIDsRef0[pidx0]
            * rtime    : the  epoch time date [s] we want to select upon! (we keep the buoy which time is closest to this `rtime`)

    RETURN:
            Updated (or not) `pidx0` and its length

            `pIDs`, `ptime` and `pidx0` have the same length!

    '''
    Nok0 = len(pIDs)
    if len(ptime)!=Nok0:
        print('ERROR [SuppressMulitOccurences]: `len(ptime)!=len(pIDs)`!'); exit(0)
    if len(pidx0)!=Nok0:
        print('ERROR [SuppressMulitOccurences]: `len(pidx0)!=len(pIDs)`!'); exit(0)
    #
    _, idxU = np.unique(pIDs, return_index=True)
    NokU = len(idxU) ; # NokU is the number of buoys once multi-occurences are removed!
    np2rm  = Nok0-NokU
    if np2rm>0:
        if iverbose>0: print('  |SMO| => ',NokU,'unique buoy IDs in an array that contains',Nok0,' buoy IDs! => '+str(np2rm)+' points to remove!')
        #
        idxOKU = pidx0[idxU] ; # because `idxU` are indices in the `pIDs` world, not in the original `pIDsRef0` world...
        # Indices of the multi-occurences:
        idxMO = np.setdiff1d( pidx0, idxOKU ) ; # keep the values of `pidx0` that are not in `idxOKU`
        del idxOKU
        zIDsMO = np.unique( pIDsRef0[idxMO] ); # unique IDs of the buoys that exist more than once in current time bin...
        #
        # Analysis:
        idxRMall = []
        for jMO in zIDsMO:
            (idxMlt,) = np.where(pIDs==jMO)
            # We keep the point that has the time the nearest to that of the center of the bin:
            jk = np.argmin(np.abs(ptime[idxMlt]-rtime))
            jKeep = idxMlt[jk]
            idxRM = np.setdiff1d( idxMlt, [jKeep] ) ; # exclude `jKeep` from idxMlt
            idxRM = pidx0[idxRM]; # translate in the `pidx0` frame!!!! IMPORTANT !!!!
            idxRMall.extend(idxRM)
        #
        idxRMall = np.array(idxRMall, dtype=int)
        if len(idxRMall)!=np2rm:
            print('ERROR [SuppressMulitOccurences]: `len(idxRMall)!=np2rm`', len(idxRMall), np2rm)
            exit(0)
        # Finally, update `pidx0`:
        pidx0 = np.setdiff1d( pidx0, idxRMall ) ; # keep the values of `pidx0` that are not in `idxRM`
        del idxRM, jKeep, idxRMall
        Nok0 = len(pidx0)
    else:
        if iverbose>0: print('    |SMO|  => no suppression to perform...')
    #
    return Nok0, pidx0




def Geo2CartNPSkm1D( pcoorG ):
    '''
         => from Geo coor. (lon,lat)[degrees] to cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj!
    '''
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    (_,n2) = np.shape(pcoorG)
    if n2!=2:
        print(' ERROR [Geo2CartNPSkm1D()]: input array `pcoorG` has a wrong a shape!')
        exit(0)
    #
    crs_src = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    crs_trg = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    #
    zx,zy,_ = crs_trg.transform_points(crs_src, pcoorG[:,1], pcoorG[:,0]).T
    #
    return np.array([ zy/1000., zx/1000. ]).T


def CartNPSkm2Geo1D( pcoorC ):
    '''
         => from cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj to Geo coor. (lon,lat)[degrees] !
    '''
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    (_,n2) = np.shape(pcoorC)
    if n2!=2:
        print(' ERROR [CartNPSkm2Geo1D()]: input array `pcoorC` has a wrong a shape!')
        exit(0)    
    #
    crs_src = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    crs_trg = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    #
    zlon,zlat,_ = crs_trg.transform_points(crs_src, 1000.*pcoorC[:,1], 1000.*pcoorC[:,0]).T
    #
    return np.array([zlat, zlon]).T




def ConvertGeo2CartesianNPSkm( plat, plon ):
    '''
         => from Geo coor. (lon,lat)[degrees] to cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj!
    '''
    #
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    crs_src = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    crs_trg = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    #
    ndim = len(np.shape(plon))
    #
    zx,zy,_ = crs_trg.transform_points(crs_src, plon, plat).T
    #
    if ndim==2:
        zx,zy = zx.T, zy.T
        #
    return zy/1000., zx/1000.



def ConvertCartesianNPSkm2Geo( pY, pX ):
    '''
         => from cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj to Geo coor. (lon,lat)[degrees] !
    '''
    #
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    crs_src = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    crs_trg = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    #
    ndim = len(np.shape(pX))
    #
    zlon,zlat,_ = crs_trg.transform_points(crs_src, 1000.*pX, 1000.*pY).T
    #
    if ndim==2:
        zlon,zlat = zlon.T, zlat.T
        #
    return zlat, zlon



def CheckTimeConsistencyQuads( kF, QD, time_dev_from_mean_allowed, iverbose=0 ):
    '''
        We look at the position DATE of all the points composing a group of quads
        (Quad class `QD`) and return a message error if some points are too far 
        in time from the mean of all the points...

        * kF: file number
        * QD: Quad class loaded from file `kF`
        * time_dev_from_mean_allowed: maximum time deviation from the mean allowed
    '''
    if iverbose>0:
        from climporn import epoch2clock
        print('\n *** In file #'+str(kF)+':')
    #
    if len(np.shape(QD.PointTime))>1:
        print('ERROR [CheckTimeConsistencyQuads()]: wrong shape for time array of Quad! (should be 1D!) => ',np.shape(QD.PointTime))
        exit(0)
    if iverbose>0: print('     (time_dev_from_mean_allowed =', time_dev_from_mean_allowed/60.,' minutes)' )
    rTmean = np.mean(QD.PointTime)
    if iverbose>0: cTmean = epoch2clock(rTmean)
    rStdDv = StdDev(rTmean, QD.PointTime)
    if iverbose>0: print('     => Actual Mean time =',cTmean,', Standard Deviation =',rStdDv/60.,' minutes!')
    zadiff = np.abs(QD.PointTime-rTmean)
    zdt = np.max(zadiff)/60.
    if iverbose>0: print('     ==> furthest point is '+str(round(zdt,2))+' minutes away from mean!')
    #
    if np.any(zadiff>time_dev_from_mean_allowed):
        print('ERROR [CheckTimeConsistencyQuads()]: some points read in file #'+str(kF)+' are too far (in time) from mean time!')
        exit(0)
    else:
        if iverbose>0: print('     => ok! No points further than '+str(time_dev_from_mean_allowed/60.)+' minutes from mean...')
    #
    return rTmean




def StreamTimeSanityCheck( cstrm, ptim, pVTb, pmsk, pBpR, iverbose=0):
    # Now, in each record of the stream we should exclude buoys which time position is not inside the expected time bin
    # or is just too far away from the mean of all buoys
    # => if such a buoy is canceld at stream # k, it should also be canceled at following records
    (NRmax,_) = np.shape(pmsk)
    zmsk = pmsk.copy()
    zBpR = pBpR.copy()
    #for jr in range(NRmax):
    #    print('SUMMARY BEFORE/ rec.',jr,': zBpR[jr], sum(zmsk[jr,:]) =',zBpR[jr], np.sum(zmsk[jr,:]))
    #
    kFU = 0
    if iverbose>0: print('\n *** Time location sanity test and fix for stream #'+str(cstrm))
    for jrec in range(NRmax):
        # At this record, the time position of all buoys of this stream is: ztim[jrec,:]
        t_mean = np.mean(ptim[jrec,:])
        rStdDv = StdDev(t_mean, ptim[jrec,:])
        zadiff = np.abs(ptim[jrec,:]-t_mean)
        zdt = np.max(zadiff)/3600.
        if iverbose>0:
            print('  * rec #',jrec,'of this stream:')            
            print('    mean time for this record is:',epoch2clock(t_mean))
            print('    bin center time, and bounds:',epoch2clock(pVTb[jrec,0]),epoch2clock(pVTb[jrec,1]),epoch2clock(pVTb[jrec,2]))
            print('    standard Deviation =',round(rStdDv/3600.,3),' hours!, nb of buoys ='+str(np.sum(zmsk[jrec,:])))
            print('     ==> furthest point is '+str(round(zdt,2))+'h away from mean! Max dev. allowed =',round(max_t_dev_allowed_in_bin/3600.,2),'h')
        #
        # Outside of the bin?
        idx_rmA = []
        lcancel = np.any( ptim[jrec,:]<pVTb[jrec,1] ) or np.any( ptim[jrec,:] > pVTb[jrec,2] )
        if lcancel:
            (idx_rm_m,) , (idx_rm_p,) = np.where( ptim[jrec,:]<pVTb[jrec,1] ) ,  np.where( ptim[jrec,:]>pVTb[jrec,2] )                
            idx_rmA = np.concatenate([ idx_rm_m , idx_rm_p ])

            if np.sum(zmsk[jrec:,idx_rmA])>0:
                #print('INSIDE before / rec.',jrec,': zBpR[jrec], sum(zmsk[jrec,:]) =',zBpR[jrec], np.sum(zmsk[jrec,:]))
                #print('INSIDE before / rec.',jrec+1,': zBpR[jrec+1], sum(zmsk[jrec+1,:]) =',zBpR[jrec+1], np.sum(zmsk[jrec+1,:]))                
                if iverbose>0:
                    print('    WARNING: the time position of some buoys are outside of that of the expected time bin!!!')
                    print('    ==> we have to cancel '+str(len(idx_rmA))+' points / '+str(np.sum(zmsk[jrec,:])))
                jr = jrec
                if jrec<2:                    
                    kFU += 1 ; # => means fields will be shrinked later on...
                    jr=0  ; # if 2nd record (jrec=1) to be canceled then the 1st record becomes useless!
                #zBpR[jr:] = zBpR[jr:] - len(idx_rmA)
                zmsk[jr:,idx_rmA] = 0 ; # This and following records!!!                
                #print('INSIDE after / rec.',jrec,': zBpR[jrec], sum(zmsk[jrec,:]) =',zBpR[jrec], np.sum(zmsk[jrec,:]))
                #print('INSIDE after / rec.',jrec+1,': zBpR[jrec+1], sum(zmsk[jrec+1,:]) =',zBpR[jrec+1], np.sum(zmsk[jrec+1,:]))

        #
        # Too far from the mean?
        idx_rmB = []
        lcancel = np.any(zadiff>max_t_dev_allowed_in_bin)
        if lcancel:
            (idx_rmB,) = np.where( zadiff>max_t_dev_allowed_in_bin )
            if np.sum(zmsk[jrec:,idx_rmB])>0:
                if iverbose>0:
                    print('    WARNING: the time position of some buoys are too far from time mean of all buoys of this bin!!!')
                    print('    ==> we have to cancel '+str(len(idx_rmB))+' points / '+str(np.sum(zmsk[jrec,:])))
                jr = jrec
                if jrec<2:                    
                    kFU += 1 ; # => means fields will be shrinked later on...
                    jr=0  ; # if 2nd record (jrec=1) to be canceled then the 1st record becomes useless!
                #zBpR[jr:] = zBpR[jr:] - len(idx_rmB)
                zmsk[jr:,idx_rmB] = 0
    #
    for jr in range(NRmax):
        zBpR[jr] = np.sum(zmsk[jr,:])
        #print('SUMMARY AFTER/ rec.',jr,': zBpR[jr], sum(zmsk[jr,:]) =',zBpR[jr], np.sum(zmsk[jr,:]))
        

    #for jr in range(NRmax):
    #    if zBpR[jr] != np.sum(zmsk[jr,:]):
    #        print('ERROR [StreamTimeSanityCheck()]: `zBpR[jr] != np.sum(zmsk[jr,:])`',zBpR[jr], np.sum(zmsk[jr,:]))
    #        exit(0)
        
    return kFU, zmsk, zBpR
                
