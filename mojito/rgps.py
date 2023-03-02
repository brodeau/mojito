import numpy as np

FillValue = -9999.

def batchSummaryRGPS( pNBini, pTcini, pIDs, pNRc ):
    from climporn import epoch2clock
    #
    (Nbtch,NbMax) = np.shape(pIDs)
    if Nbtch != len(pNBini):
        print('ERROR: [batchSummaryRGPS()] => error #1')
    if not NbMax == max(pNBini):
        print('ERROR: [batchSummaryRGPS()] => error #2')
    print('\n ==========   SUMMARY   ==========')
    print(' *** Number of identified batches: '+str(Nbtch))
    print(' *** Number of buoys selected in each batch:')
    for js in range(Nbtch):
        cTc0  = epoch2clock(pTcini[js])
        print('        * Batch #'+str(js)+' initiated at time bin centered around '+cTc0+' => has '+str(pNBini[js])+' buoys')
    print(' *** Max number of buoys possibly found in a batch = ',NbMax)
    print('     * shape of ZIDs =', np.shape(pIDs))
    print('     * shape of ZNRc =', np.shape(pNRc))
    print(' ===================================\n')


def KeepDataInterest( dt1, dt2, ptime0, pIDs0 ):
    '''
       Only keep data of interest in 1D arrays (except time array), based on date1 and date2.
       Excluded points are masked...
       Again, time array `ptime0` won't be masked (unnecessary and dangerous)
       Input:
                * dt1, dt2 : start & end time ([s] UNIX epoch time)

       Returns:
                * nB   : number of unique buoys that exist for at least 1 record during specified date range
                * zIDs : list (unique) of IDs for these buoys (array[nB] of int)
                * idx_msk: indices to mask to retain only data of interest (date range) in raw data arrays
    '''
    (nP,) = np.shape( pIDs0 )
    zIDs0 = pIDs0.copy()
    #
    # Will mask all point that are before and beyond our period of interest:
    (idx_msk,) = np.where( (ptime0 < dt1) | (ptime0 > dt2) )
    (idx_fub,) = np.where(  zIDs0 < 0 ) ; # just in case some negative IDs exist
    idx_msk    = np.unique( np.concatenate([idx_msk,idx_fub]) )
    #
    zIDs0[idx_msk] = -999
    #
    # Remaining buoys (IDs)
    (idx,) = np.where(zIDs0 >= 0) ; # yes, ID=0 exists in RGPS data
    zIDs = np.sort( np.unique( zIDs0[idx] ) ) ; # if not `[idx]` then `rmskVal` is counted once!!!
    nB   = len(zIDs)
    print('   * [KDI] We found '+str(nB)+' unique buoys alive during specified period of time.')
    #
    return nB, zIDs, idx_msk



def LoadData4TimeRange( idate1, idate2, fRGPS, l_doYX=False ):
    '''

    RETURNS:
     * nPr: number of points of interst
     * nBu: number of unique buoys of interest
     * zIDsU: unique IDs of the buoys of interest
    
    '''
    from .ncio import LoadDataRGPS

    # Open, inspect the input file and load raw data:
    if l_doYX:
        Np0, _, ztime0, zykm0, zxkm0, zlat0, zlon0, zIDs0, _ = LoadDataRGPS( fRGPS )
    else:
        Np0, _, ztime0,     _,     _, zlat0, zlon0, zIDs0, _ = LoadDataRGPS( fRGPS )
    print('   * [LD4TR] before time-range-exclusion we have '+str(Np0)+' points in the file.')

    nBu, zIDsU, idxmsk = KeepDataInterest( idate1, idate2, ztime0, zIDs0 )
    # * nBu: number of unique buoys that exist for at least 1 record during specified date range aka whole period (WP)
    # * zIDsU : array(nBu) list (unique) of IDs for these buoys
    # * idxmsk: indices of points to cancel

    # Masking all point that are before and beyond our period of interest (note: `ztime0` is not masked!):
    zmsk0 = np.ones(Np0, dtype='i1')
    zmsk0[idxmsk] = 0
    lcond = ( zmsk0 == 0 )
    zIDs0[idxmsk] = int(FillValue) ; zIDs0 =  np.ma.masked_where( lcond, zIDs0 )
    zlat0[idxmsk] =     FillValue  ; zlat0 =  np.ma.masked_where( lcond, zlat0 )
    zlon0[idxmsk] =     FillValue  ; zlon0 =  np.ma.masked_where( lcond, zlon0 )
    if l_doYX:
        zykm0[idxmsk] = FillValue  ; zykm0 =  np.ma.masked_where( lcond, zykm0 )
        zxkm0[idxmsk] = FillValue  ; zxkm0 =  np.ma.masked_where( lcond, zxkm0 )
    nPr = np.sum(zmsk0)
    print('   * [LD4TR] after time-range-exclusion we have '+str(nPr)+' / '+str(Np0)+' points.')
    del zmsk0, lcond
    
    if l_doYX:
        return nPr, nBu, zIDsU, ztime0, zIDs0, zlat0, zlon0, zykm0, zxkm0
    else:
        return nPr, nBu, zIDsU, ztime0, zIDs0, zlat0, zlon0


def ValidNextRecord( time_min, kidx, ptime0, pBIDs0, pidxIgnore, dtNom, devdtNom ):
    '''
       RETURNS:
            * nbROK: number of ok records for this point
                     0 -> this buoy (kID) has a mono record in the whole period (not only in the bin) => should be canceled
                     1 -> this point (kidx) has successors in time but none is reasonably timed (based on dtNom & devdtNom)
                     2 -> a reasonably timed successor has been found
    '''
    from climporn import epoch2clock
    #
    idxScsr = -9999
    zt0     = -9999.
    zts     = -9999.
    nbROK   = 0
    #print( 'LOLO: [ValidNextRecord] 0: dtNom, devdtNom =',dtNom/3600, devdtNom/3600)

    kID = pBIDs0[kidx]
    (idxBuoy,) = np.where( pBIDs0 == kID )
    #print( 'LOLO: [ValidNextRecord] 1: we have '+str(len(idxBuoy))+' buoy positions,',str(len(np.unique(idxBuoy))) )
    #print(' LOLO: dates for all these positions:')
    #print(' Buoy ID =',kID)
    #print(' idxBuoy =',idxBuoy,' time:')
    #for ki in idxBuoy:
    #    print(epoch2clock(ptime0[ki]),end=' ')
    #print('')
    
    if len(idxBuoy) > 1:
        # Ok, now we know that this point is not the mono-occurence of a mono-record buoy (in the whole period not only in the bin)
        # We focus on time location of this buoy after the begining of current bin:
        (idxBuoy,) = np.where( (ptime0>=time_min) & (pBIDs0 == kID) )
        #
        zt0  = ptime0[idxBuoy[0]] ; # the first time position in this bin
        ztR1 = zt0 + dtNom - devdtNom ; # Reasonable lower time bond for successor point
        ztR2 = zt0 + dtNom + devdtNom ; # Reasonable lower time bond for successor point
        nbROK   = 1
        lHasSuccessor = np.any( (ptime0[idxBuoy]>ztR1) & (ptime0[idxBuoy]<ztR2) )
        #
        if lHasSuccessor:
            nbROK   = 2
            (idxS,) = np.where( (ptime0[idxBuoy]>ztR1) & (ptime0[idxBuoy]<ztR2) )
            if len(idxS)!=1:
                #print(' [ValidNextRecord()]: point has more than 1 possible successor!!!')
                #print(' Point is:',epoch2clock(zt0))
                #print(' successor are:')
                #for zt in ptime0[idxBuoy[idxS]]: print(epoch2clock(zt))
                ii = np.argmin( np.abs(ptime0[idxBuoy[idxS]] - (zt0 + dtNom)) )
                idxScsr = idxS[ii]
                #print(' => selected:',epoch2clock(ptime0[idxBuoy[idxScsr]]) ); exit(0)
            else:
                idxScsr = idxS[0]
            #
            idxScsr = idxBuoy[idxScsr] ; # in the ref0 frame!
            zts = ptime0[idxScsr]

            if idxScsr in pidxIgnore:
                print('ERROR [ValidNextRecord()] `idxScsr` is forbidden by `pidxIgnore`! ')
                exit(0)            
    #
    return nbROK, np.array([kidx, idxScsr], dtype=int), np.array([zt0, zts])




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

        nP, zidx = ExcludeMultiOccurences( vids, vtime, vids_ref, vidx, t_ref )

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


def EMO2( pIDs, ptime, pIDsRef0, ptimeRef0, pidx0, binTctr, criterion='nearest', dtNom=3600*24*3, devdtNom=3600*6, iverbose=0 ):
    '''
         For many possible reasons the same buoy ID can have multiple time-position occurences in a given time bin
         (especially if wide bin)
         => we need to keep only one of these occurences of the most approriate location within this time bin
         => one of the selection process is too look which of these multi-occuring positions is more promising
            in terms of upcomming position after about 3 days (`dtNom`) !
    INPUT:
            * pIDs     : 1D array of integers containing buoy IDs with possible multi-occurence of certain IDs
            * ptime    : 1D array of real containing epoch time date associated to each buoy [s] (UNIX time)
            * pIDsRef0 : 1D array of integers containing buoy IDs the "0" reference
            * pidx0    : 1D array of integers containing indices that do this: pIDs == pIDsRef0[pidx0]
            * binTctr  : time at center of time bin/range we are dealing with [s] (UNIX time)
            * criterion: criterion used to decide which of the multi-occurence of the same buoy within the time bin/range should stay!
                         - 'nearest' => keep the occurence closest to that of center of time bin
                         - 'first'  => keep the occurence closest to begining of time bin
                         - 'last'   => keep the occurence closest to the end of time bin
                         - 'successors' => keep the occurence that is the most likely to have an appropriate successor in time

    RETURN:
            Updated (or not) `pidx0` and its length

            `pIDs`, `ptime` and `pidx0` have the same length!

    '''
    if not criterion in ['nearest','first','last']:
        print('ERROR [EMO2]: criterion "'+criterion+'" is unknown!')
        exit(0)
    #
    if np.shape(ptimeRef0) != np.shape(pIDsRef0):
         print('ERROR [EMO2]: `ptimeRef0` has a wrong shape!', np.shape(ptimeRef0))
         exit(0)
    if np.sum(np.abs(ptimeRef0[pidx0] - ptime)) != 0:
        print('ERROR [EMO2]: `ptimeRef0[pidx0]` not equal to `ptime` !')
        exit(0)
    #
    (Nok0,) = np.shape(pIDs)
    if np.shape(ptime)!=(Nok0,):
        print('ERROR [EMO2]: `np.shape(ptime)!=np.shape(pIDs)`!'); exit(0)
    if np.shape(pidx0)!=(Nok0,):
        print('ERROR [EMO2]: `np.shape(pidx0)!=np.shape(pIDs)`!'); exit(0)

    # Unique buoys ?
    _, idxU = np.unique(pIDs, return_index=True)
    NokU    = len(idxU) ; # NokU is the number of buoys once multi-occurences are removed!
    np2rm   = Nok0-NokU    
    
    if np2rm>0:
        if iverbose>0: print('    * [EMO2] => ',NokU,'unique buoy IDs in array featuring ',Nok0,' buoy IDs! => '+str(np2rm)+' points to exclude!')
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
            #
            # ----------------
            idxMlt0 = pidx0[idxMlt]
            isuccess = np.zeros(len(idxMlt),dtype='i1') ; # "1" if an acceptable successor was found for this point 0 otherwize
            #
            jm = 0
            for idx0 in idxMlt0:
                zt0 = ptimeRef0[idx0]
                ztf_ideal = zt0 + dtNom
                (idxFuture0,) = np.where( (pIDsRef0==jMO) & (ptimeRef0>zt0) )
                if len(idxFuture0)>0:
                    zdevFideal = np.abs(ptimeRef0[idxFuture0] - ztf_ideal)
                    kK = np.argmin(zdevFideal)
                    if zdevFideal[kK] < devdtNom:
                        isuccess[jm] = 1
                jm += 1
            #
            NpS = np.sum(isuccess)
            if NpS == 0:
                # No candidate, none of the multi positions of this buoy has a possible good successor in time :(
                # => keep the time position closest to the center of bin...
                # In the end it does not really matter which one we keep since no good successor
                # means it won't be used later...
                #jk = np.argmin( ptime[idxMlt] )
                jk = np.argmin(np.abs(ptime[idxMlt]-binTctr))
                #
            else:
                # Nps > 0 !
                (iis,) = np.where(isuccess==1)
                if NpS == 1:
                    # Only 1 of the multi positions of this buoy has a possible good successor in time :)
                    # => obviously the one we pick
                    jk = iis[0]
                    #
                else:
                    # More than 1 of the multi positions of this buoy has a possible good successor in time :)
                    # => we have to pick one of them based on the specified criterion
                    if   criterion=='nearest':
                        jk = iis[ np.argmin(np.abs(ptime[idxMlt[iis]]-binTctr)) ]
                    elif criterion=='first':
                        jk = iis[0]
                    elif criterion=='last':                    
                        jk = iis[-1]
            # ----------------
            #
            idxRM = np.setdiff1d( idxMlt, [idxMlt[jk]] ) ; # exclude `idxMlt[jk]` from `idxMlt`
            idxRMall.extend( pidx0[idxRM] ) ; # add to exclusion list, once converted to `pidx0` frame!
        #
        if len(idxRMall)!=np2rm:
            print('ERROR [EMO2]: `len(idxRMall)!=np2rm`', len(idxRMall), np2rm) ; exit(0)
        #
        pidx0 = np.setdiff1d( pidx0, idxRMall ) ; # exclude `idxRMall` from `pidx0`
        Nok0 = len(pidx0)
        if iverbose>0: print('    * [EMO2] => excluded '+str(np2rm)+' pt. of time bin due to multi-occur. of same buoy ID! (criterion='+criterion+')')
    else:
        if iverbose>0: print('    * [EMO2] => no suppression to perform...')
    #
    return Nok0, pidx0






















def ExcludeMultiOccurences( pIDs, ptime, pIDsRef0, pidx0, binTctr, criterion='nearest',
                            ptimeRef0=[], dtNom=None, devdtNom=None, iverbose=0 ):
    '''
         For many possible reasons the same buoy ID can exist more than once in `pIDs`,
         => we need to keep only one occurence of the location index of these points,
            based on the date (closest to center of bin `rTc`)
    INPUT:
            * pIDs     : 1D array of integers containing buoy IDs with possible multi-occurence of certain IDs
            * ptime    : 1D array of real containing epoch time date associated to each buoy [s] (UNIX time)
            * pIDsRef0 : 1D array of integers containing buoy IDs the "0" reference
            * pidx0    : 1D array of integers containing indices that do this: pIDs == pIDsRef0[pidx0]
            * binTctr  : time at center of time bin/range we are dealing with [s] (UNIX time)
            * criterion: criterion used to decide which of the multi-occurence of the same buoy within the time bin/range should stay!
                         - 'nearest' => keep the occurence closest to that of center of time bin
                         - 'first'  => keep the occurence closest to begining of time bin
                         - 'last'   => keep the occurence closest to the end of time bin
                         - 'successors' => keep the occurence that is the most likely to have an appropriate successor in time

    RETURN:
            Updated (or not) `pidx0` and its length

            `pIDs`, `ptime` and `pidx0` have the same length!

    '''
    if not criterion in ['nearest','first','last','successors']:
        print('ERROR [ExcludeMultiOccurences]: criterion "'+criterion+'" is unknown!')
        exit(0)
    #
    if criterion=='successors' and (not np.shape(ptimeRef0)==np.shape(pIDsRef0) or not dtNom or not devdtNom):
         print('ERROR [ExcludeMultiOccurences]: with "'+criterion+'" criterion, you must provide `ptimeRef0` array and `dtNom` and `devdtNom`!')
         exit(0)        
    #
    (Nok0,) = np.shape(pIDs)
    if np.shape(ptime)!=(Nok0,):
        print('ERROR [ExcludeMultiOccurences]: `np.shape(ptime)!=np.shape(pIDs)`!'); exit(0)
    if np.shape(pidx0)!=(Nok0,):
        print('ERROR [ExcludeMultiOccurences]: `np.shape(pidx0)!=np.shape(pIDs)`!'); exit(0)
    #
    _, idxU = np.unique(pIDs, return_index=True)
    NokU    = len(idxU) ; # NokU is the number of buoys once multi-occurences are removed!
    np2rm   = Nok0-NokU
    if np2rm>0:
        if iverbose>0: print('    * [EMO] => ',NokU,'unique buoy IDs in array featuring ',Nok0,' buoy IDs! => '+str(np2rm)+' points to exclude!')
        #
        idxOKU = pidx0[idxU] ; # because `idxU` are indices in the `pIDs` world, not in the original `pIDsRef0` world...
        # Indices of the multi-occurences:
        idxMO = np.setdiff1d( pidx0, idxOKU ) ; # keep the values of `pidx0` that are not in `idxOKU`
        del idxOKU
        zIDsMO = np.unique( pIDsRef0[idxMO] ); # unique IDs of the buoys that exist more than once in current time bin...

        if criterion=='successors':
            #from climporn import epoch2clock
            # Checking that `ptimeRef0` is consistent with `ptime`
            if np.sum(np.abs(ptimeRef0[pidx0] - ptime)) != 0:
                print('ERROR [ExcludeMultiOccurences]: `ptimeRef0[pidx0]` not equal to `ptime` !')
                exit(0)
        
        # Analysis:
        idxRMall = []
        for jMO in zIDsMO:
            (idxMlt,) = np.where(pIDs==jMO)
            #
            if   criterion=='nearest':
                # We keep the point that has the time the nearest to that of the center of the bin:
                jk = np.argmin(np.abs(ptime[idxMlt]-binTctr))
            elif criterion=='first':
                jk = np.argmin( ptime[idxMlt] )
            elif criterion=='last':
                jk = np.argmax( ptime[idxMlt] )
            #
            elif criterion=='successors':
                idxMlt0 = pidx0[idxMlt]
                #print(' Times for all the occurences of buoy inside bin:')
                #for zt in ptimeRef0[idxMlt0]: print(epoch2clock(zt))
                isuccess = np.zeros(len(idxMlt),dtype='i1')          ; # "1" if an acceptable successor was found for this point 0 otherwize
                #
                jm = 0
                for idx0 in idxMlt0:
                    #print(' * idx0 =', idx0,':')
                    zt0 = ptimeRef0[idx0]
                    ztf_ideal = zt0 + dtNom
                    #print('   => time =',epoch2clock(zt0), '  ==> successor ideal time =',epoch2clock(ztf_ideal))
                    (idxFuture0,) = np.where( (pIDsRef0==jMO) & (ptimeRef0>zt0) )
                    if len(idxFuture0)>0:
                        #print('   ==> future occurence of the buoy: ')
                        #for ztf in ptimeRef0[idxFuture0]: print('         * '+epoch2clock(ztf))
                        zdevFideal = np.abs(ptimeRef0[idxFuture0] - ztf_ideal)
                        #print('   ==> deviation from ideal in hours:',zdevFideal[:]/3600)
                        kK = np.argmin(zdevFideal)
                        #print( '   ==> ',zdevFideal[kK]/3600,'wins! => ', epoch2clock(ptimeRef0[idxFuture0[kK]]) )
                        if zdevFideal[kK] < devdtNom:
                            isuccess[jm] = 1
                    #
                    jm += 1
                #print(' isuccess =',isuccess)
                if np.sum(isuccess) == 0:
                    # No candidate, keeping earliest point, just as in 'first' criterion
                    jk = np.argmin( ptime[idxMlt] )
                else:
                    (iis,) = np.where(isuccess==1)
                    jk = iis[0] ; # keeping first of the winners...                    
            #   
            jKeep = idxMlt[jk]
            idxRM = np.setdiff1d( idxMlt, [jKeep] ) ; # exclude `jKeep` from idxMlt
            idxRM0 = pidx0[idxRM]; # translate in the `pidx0` frame!!!! IMPORTANT !!!!
            idxRMall.extend(idxRM0)
            #
            #print('FINAL we keep time:', epoch2clock(ptimeRef0[pidx0[jKeep]]),'!')
            #for idx in idxRM0:
            #    print('FINAL we delete time:', epoch2clock(ptimeRef0[idx]),'!\n')            
        #
        idxRMall = np.array(idxRMall, dtype=int)
        if len(idxRMall)!=np2rm:
            print('ERROR [ExcludeMultiOccurences]: `len(idxRMall)!=np2rm`', len(idxRMall), np2rm)
            exit(0)
        # Finally, update `pidx0`:
        pidx0 = np.setdiff1d( pidx0, idxRMall )
        Nok0 = len(pidx0)
        if iverbose>0: print('    * [EMO] => excluded '+str(np2rm)+' pt. of time bin due to multi-occur. of same buoy ID! (criterion='+criterion+')')
    else:
        if iverbose>0: print('    * [EMO]  => no suppression to perform...')
    #
    return Nok0, pidx0



def BatchTimeSanityCheck( cbtch, ptim, pVTb, pmsk, pBpR, tdev_max, iverbose=0):
    from .util import StdDev
    '''
       * cbtch: string to identufy current batch
       * ptim:  2D masked time array  (NCRmax,NvB) or original RGPS time positions [s]
       * pVTb:  "almost 1D" time array (NCRmax,3) of time bins used [s]
       * pmsk:  2D mask array (NCRmax,NvB) for canceled points (int)
       * pBpR:  1D (NCRmax) array of the number of buoys alive at each record (int)
       * tdev_max: max authorized time deviation from the mean for a buoy [s]
    '''

    
    # Now, in each record of the batch we should exclude buoys which time position is not inside the expected time bin
    # or is just too far away from the mean of all buoys
    # => if such a buoy is canceld at batch # k, it should also be canceled at following records
    (NRmax,_) = np.shape(pmsk)
    zmsk = pmsk.copy()
    zBpR = pBpR.copy()
    #for jr in range(NRmax):
    #    print('SUMMARY BEFORE/ rec.',jr,': zBpR[jr], sum(zmsk[jr,:]) =',zBpR[jr], np.sum(zmsk[jr,:]))
    #
    kFU = 0
    if iverbose>0: print('  * [BatchTimeSanityCheck]: time location sanity test and fix for batch #'+str(cbtch))
    for jrec in range(NRmax):
        # At this record, the time position of all buoys of this batch is: ztim[jrec,:]
        t_mean = np.mean(ptim[jrec,:])
        rStdDv = StdDev(t_mean, ptim[jrec,:])
        zadiff = np.abs(ptim[jrec,:]-t_mean)
        zdt = np.max(zadiff)/3600.
        if iverbose>0:
            from climporn import epoch2clock
            print('  * rec #',jrec,'of this batch:')
            print('    mean time for this record is:',epoch2clock(t_mean))
            print('    bin center time, and bounds:',epoch2clock(pVTb[jrec,0]),epoch2clock(pVTb[jrec,1]),epoch2clock(pVTb[jrec,2]))
            print('    standard Deviation =',round(rStdDv/3600.,3),' hours!, nb of buoys ='+str(np.sum(zmsk[jrec,:])))
            print('     ==> furthest point is '+str(round(zdt,2))+'h away from mean! Max dev. allowed =',round(tdev_max/3600.,2),'h')
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
        lcancel = np.any(zadiff>tdev_max)
        if lcancel:
            (idx_rmB,) = np.where( zadiff>tdev_max )
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
    #        print('ERROR [BatchTimeSanityCheck()]: `zBpR[jr] != np.sum(zmsk[jr,:])`',zBpR[jr], np.sum(zmsk[jr,:]))
    #        exit(0)

    return kFU, zmsk, zBpR



