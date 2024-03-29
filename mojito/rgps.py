import numpy as np

dt0_RGPS = 3*24*3600 ; # the expected nominal time step of the input data, ~ 3 days [s]

from .ncio import FillValue


def batchSummaryRGPS( pTbin, pNBini, pidxBinN, pIDs, pNRc ):
    from .util import epoch2clock as ep2c
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
        cTc0, cTcA, cTcB  = ep2c(pTbin[pidxBinN[js],0]), ep2c(pTbin[pidxBinN[js],1]), ep2c(pTbin[pidxBinN[js],2])
        print('        * Batch #'+str(js)+' initiated at time bin centered around '+cTc0+' => has '+str(pNBini[js])+' buoys')
        print('                         => bounds of time bin: '+cTcA+' -- '+cTcB)
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


def ValidUpComingRecord( time_min, kidx, ptime0, pBIDs0, devdtNom=6*3600 ):
    '''
        * time_min: date at the lower time bound of the considered time bin [s]
        * kidx:     integer => the index that accesses the currently treated buoy in the `*0` original arrays like `ptime0` or `pBIDs0`!
        * ptime0:   unmasked original time array data
        * pBIDs0:   masked original array of IDs (it is masked outside the whole period of interest)

       RETURNS:
            * nbROK: number of ok records for this point
                     0 -> this buoy (kID) has a mono record in the whole period (not only in the bin) => should be canceled
                     1 -> this point (kidx) has successors in time but none is reasonably timed (based on dt0_RGPS & devdtNom)
                     2 -> a "reasonably-timed" upcoming position has been found for this buoy
            * zidx0_recs: valid records including this one aka `kidx` in first position

        "VUCR" => Valid UpComing Record !
    '''
    from .util import epoch2clock as ep2c ; #debug
    #
    idxUC = -9999
    zt0   = -9999.
    nbROK = 0
    zidx0_recs = []

    kID = pBIDs0[kidx]                     ; # the ID of the buoy we are dealing with

    #(idxBuoy,) = np.where( (ptime0>=time_min) & (pBIDs0 == kID) ) ; # => probably more demanding than bellow:
    (idxBuoy,) = np.where( pBIDs0 == kID ) ;              # restricts to the buoy of interest
    (idxT,)    = np.where( ptime0[idxBuoy] >= time_min ); # restricts to the time of interest
    idxBuoy    = idxBuoy[idxT]
    del idxT

    # If time bins are wide, `idxBuoy` can contain indices before (smaller than) `kidx`, so clean this:
    if len(idxBuoy) > 0:
        if idxBuoy[0]<kidx:
            (idxKeep,) = np.where(idxBuoy>=kidx)
            idxBuoy = idxBuoy[idxKeep]
            del idxKeep

    if len(idxBuoy) > 1:
        # Ok, now we know that this point is not the mono-occurence of a mono-record buoy (in the whole period not only in the bin)
        nbROK = 1
        #
        idxT0 = idxBuoy[0]
        zt0   = ptime0[idxT0] ; # the first time position inside this bin
        ztUC  = zt0 + dt0_RGPS
        ztR1  = ztUC - devdtNom ; # Reasonable lower time bond for successor point
        ztR2  = ztUC + devdtNom ; # Reasonable upper time bond for successor point
        #
        ll1 = (ptime0[idxBuoy]>ztR1)
        ll2 = (ptime0[idxBuoy]<ztR2)
        lHasVUCR = np.any( ll1 & ll2 ) ; # 1st record is ignored by definition
        #
        if lHasVUCR:
            # There is at least 1 upcoming position of our buoy within the "reasonable" expected upcomming time bin
            nbROK   = 2
            (idxS,) = np.where( ll1 & ll2 )
            if len(idxS)>1:
                #print('LOLOLOLOLOLO: more than one buoy position found in upcoming time bin!!! n=',len(idxS))
                #print('LOLO: time of buoy:',ep2c(zt0))
                #print('LOLO: expected ideal time is:',ep2c(ztUC))
                #print('LOLO: these times are:')
                #for it in ptime0[idxBuoy[idxS]]: print('       =>',ep2c(it))                
                # There are actually more than 1! Must pick the closest to expected time:
                ii = np.argmin( np.abs(ptime0[idxBuoy[idxS]] - ztUC) )
                idxUC = idxS[ii]
                #print('  ==> selected:',ep2c(ptime0[idxBuoy[idxUC]]))
                #print('')
            else:
                idxUC = idxS[0]
            #
            idxUC = idxBuoy[idxUC] ; # in the ref0 frame!

            zidx0_recs = np.array([idxT0, idxUC], dtype=int)
        else:
            zidx0_recs = np.array([idxT0], dtype=int)

        if zidx0_recs[0] != kidx:
            print('ERROR: [ValidUpComingRecord] => something wrong! `zidx0_recs[0] != kidx`')
            print('  => kidx       =', kidx)
            print('  => idxBuoy    =', idxBuoy)
            print('  => zidx0_recs =',zidx0_recs)
            exit(0)

    ### if len(idxBuoy) > 1
    #
    return nbROK, zidx0_recs




def ValidCnsctvRecordsBuoy( time_min, kidx, ptime0, pBIDs0, pidx_ignore, dt_expected, devdtNom ):
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
                * devdtNom: overshoot tolerance for `dt_expected` [s]
       Returns:
               * nbROK   : number of valid consecutive records for this buoy
               * idx0_id : array of location indices (in the raw data arrays) for these valid records of this buoy
               * ztime   : array of dates associated with all these records [s]
    '''
    from .util import epoch2clock as ep2c
    #
    (idx_buoy,   ) = np.where( pBIDs0 == pBIDs0[kidx] )
    (idx_exclude,) = np.where( ptime0 < time_min )
    idx_keep0      = np.setdiff1d( idx_buoy,  idx_exclude ) ; # keep the values of `idx_buoy` that are not in `idx_exclude`
    idx_keep       = np.setdiff1d( idx_keep0, pidx_ignore ) ; # keep the values of `idx_keep0` that are not in `pidx_ignore`
    del idx_exclude, idx_buoy, idx_keep0
    #
    ztime  = ptime0[idx_keep] ; # all time records for this particular buoy
    if np.any( ztime<=time_min ):
        print('ERROR: [ValidCnsctvRecordsBuoy] => `np.any( ztime<=time_min )` !!!!'); exit(0)
    #
    nbR1b = len(idx_keep)      ; # n. of time records for this particulat buoy
    #
    nbROK = nbR1b
    #
    # We cut the series at the first occurence of a dt (time between 2 consec. pos.) > `devdtNom`
    ztime_ideal = np.array( [ ztime[0]+float(i)*float(dt_expected) for i in range(nbR1b) ], dtype=float )
    vtdev = np.abs(ztime - ztime_ideal)
    if np.any(vtdev > devdtNom):
        #print('LOLO: vtdev[:]/3600. = ',vtdev[:]/3600.)
        (indFU,) = np.where(vtdev > devdtNom)
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
    from .util import epoch2clock as ep2c

    nbf = len(list_npz_files)

    # First round to check the number of points in each file
    vit, vNbP = [], []
    for cf in list_npz_files:
        with np.load(cf) as data:
            it = int( data['itime'] ) ; # because otherwize it's an array of size 1 ?
            vit.append(it)
            npf = len(data['vids'])
            vNbP.append(npf)
        if iverbose>0: print('  '+cf+' => time ='+ep2c(it)+' | '+str(npf)+' points!')
    #
    # For merged file:
    itime_mean = int( round(np.mean(vit),0) )
    nP = np.sum(vNbP)
    if iverbose>0: print(' vNbP =', vNbP,'=>',nP,'points in total!, time =',ep2c(itime_mean))

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

    if ep2c(itime_mean) != ep2c(t_ref):
        print(' WARNING: [util.mergeNPZ] => `epoch2clock(itime_mean) != epoch2clock(t_ref)` !!!')

    # Time to save in the new npz file:
    if iverbose>0: print('  [util.mergeNPZ] ==> saving merged files into '+cf_out+' !')
    np.savez_compressed( cf_out, itime=itime_mean, date=ep2c(itime_mean), Npoints=nP, vids=vids,
                         vtime=vtime, vx=vx, vy=vy, vlon=vlon, vlat=vlat,  )
    return 0





def EMO2( pIDs, ptime, pIDsRef0, ptimeRef0, pidx0, binTctr, criterion='nearest', devdtNom=3600*6, iverbose=0 ):
    '''

                        "  Exclude Muti-Occurences of a given buoy under a given time period"

         For many possible reasons the same buoy ID can have multiple time-position occurences in a given time bin
         (especially in the case of a long time bin, > 3 days)
         => we need to keep only one of these occurences of the most approriate location within this time bin
         => one of the selection process is too look which of these multi-occuring positions is more promising
            in terms of upcomming position after about 3 days (`dt0_RGPS`) !
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
                ztf_ideal = zt0 + dt0_RGPS
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






def BatchTimeSanityCheck( cbtch, ptim, pmsk, pBpR, tdev_max, pTbin_bounds,  iverbose=0):
    from .util import StdDev
    '''
       * cbtch: string to identufy current batch
       * ptim:  2D masked time array  (NRmax,NvB) or original RGPS time positions [s]
       * pmsk:  2D mask array (NRmax,NvB) for canceled points (int)
       * pBpR:  1D (NRmax) array of the number of buoys alive at each record (int)
       * tdev_max: max authorized time deviation from the mean for a buoy [s]
       * pTbin_bounds: array (NRmax,2) of bounds of time bin in which each record was IDed (2 records can be in the same time bin if dt_bin was big !!!)
    '''
    # Now, in each record of the batch we should exclude buoys which time position is not inside the expected time bin
    # or is just too far away from the mean of all buoys
    # => if such a buoy is canceld at batch # k, it should also be canceled at following records
    (NRmax,nB) = np.shape(ptim)
    #if np.shape(pVTb)!=(NRmax,3):
    #    print('ERROR [BatchTimeSanityCheck]: wrong shape for `pVTb`',np.shape(pVTb)); exit(0)
    if np.shape(pmsk)!=(NRmax,nB):
        print('ERROR [BatchTimeSanityCheck]: wrong shape for `pmsk`',np.shape(pmsk)); exit(0)
    if np.shape(pBpR)!=(NRmax,):
        print('ERROR [BatchTimeSanityCheck]: wrong shape for `pBpR`',np.shape(pBpR)); exit(0)
    #
    zmsk = pmsk.copy()
    zBpR = pBpR.copy()
    #for jr in range(NRmax):
    #    print('SUMMARY BEFORE/ rec.',jr,': zBpR[jr], sum(zmsk[jr,:]) =',zBpR[jr], np.sum(zmsk[jr,:]))
    #
    kFU = 0
    if iverbose>0: print('  * [BTSC]: time location sanity test and fix for batch #'+str(cbtch))
    for jrec in range(NRmax):
        # At this record, the time position of all buoys of this batch is: ztim[jrec,:]
        t_mean = np.mean(ptim[jrec,:])
        rStdDv = StdDev(t_mean, ptim[jrec,:])
        zadiff = np.abs(ptim[jrec,:]-t_mean)
        zdtworse = np.max(zadiff)/3600.
        if iverbose>0:
            from .util import epoch2clock as ep2c
            print('  * [BTSC]: rec #',jrec,'of this batch:')
            print('  * [BTSC]: mean time for this record is:',ep2c(t_mean))
            print('  * [BTSC]: bounds of time bin used:',ep2c(pTbin_bounds[jrec,0]), ep2c(pTbin_bounds[jrec,1]))
            #print('  * [BTSC]: bin used (lb,c,hb) =>',ep2c(pVTb[jrec,1])+' | '+ep2c(pVTb[jrec,0])+' | '+ep2c(pVTb[jrec,2]))
            print('  * [BTSC]: standard Deviation =',round(rStdDv/3600.,3),' hours!, nb of buoys ='+str(np.sum(zmsk[jrec,:])))
            print('  * [BTSC]:  ==> furthest point is '+str(round(zdtworse,2))+'h away from mean! Max dev. allowed =',round(tdev_max/3600.,2),'h')
        
        ## Outside of the bin?
        if t_mean<pTbin_bounds[jrec,0] or t_mean>pTbin_bounds[jrec,1]:
            print(' ERROR [BatchTimeSanityCheck]: mean time for selected points does not fall inside expected bin!!!')
            exit(0)
        ## => mind that with wide time bins the second record can still be inside the first bin !!!
        ##zoutside = (ptim[jrec,:]<pVTb[jrec,1]) | (ptim[jrec,:] > pVTb[jrec,2])
        #zoutside = (ptim[jrec,:]<pVTb[0,1]) | (ptim[jrec,:] > pVTb[jrec,2]) ; # => using 0 rather than jrec in `pVTb` !!!
        #lOutside = np.any( zoutside )
        #if lOutside:
        #    print(' ERROR [BatchTimeSanityCheck]: the time position of some buoys are outside of what seems reasonable!!!')
        #    (idx_rmO,) = np.where( zoutside )
        #    print('  * [BTSC]: =>  jrec, ptim[jrec,idx_rmO] =', jrec, np.array( [ ep2c(ptim[jrec,i]) for i in idx_rmO ] ) )
        #    print('  * [BTSC]: => pVTb[0,1], pVTb[jrec,1]',ep2c(pVTb[0,1]), ep2c(pVTb[jrec,1]))
        #    print('  * [BTSC]: =>     pVTb[jrec,2]', ep2c(pVTb[jrec,2]))
        #    exit(0)
        #    #(idx_rmO,) = np.where( zoutside )
        #    #if np.sum(zmsk[jrec:,idx_rmO])>0:
        #    #    if iverbose>0:
        #    #        print(' WARNING [BatchTimeSanityCheck]: the time position of some buoys are outside of that of the expected time bin!!!')
        #    #        print('  * [BTSC]: ==> we have to cancel '+str(len(idx_rmO))+' points / '+str(np.sum(zmsk[jrec,:])))
        #    #    jr = jrec
        #    #    if jrec<2:
        #    #        kFU += 1 ; # => means fields will be shrinked later on...
        #    #        jr=0  ; # if 2nd record (jrec=1) to be canceled then the 1st record becomes useless!
        #    #    zmsk[jr:,idx_rmO] = 0 ; # This and following records!!!

        # Too far from the mean?
        idx_rmB = []
        lTooFar = np.any(zadiff>tdev_max)
        if lTooFar:
            (idx_rmB,) = np.where( zadiff>tdev_max )
            if np.sum(zmsk[jrec:,idx_rmB])>0:
                if iverbose>0:
                    print('  * [BTSC]: WARNING: the time pos. of some buoys too far from time mean of all buoys of this record of this batch!')
                    print('  * [BTSC]: ==> we have to cancel '+str(len(idx_rmB))+' points / '+str(np.sum(zmsk[jrec,:])))
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

