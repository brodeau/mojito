#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J xHORRORx
#SBATCH -o out_RGPS_selection_%J.out
#SBATCH -e err_RGPS_selection_%J.err
#SBATCH --time=02:55:00
#SBATCH --account=python
#SBATCH --mem=8000
##################################################################

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np
from re import split

import mojito as mjt

idebug = 0
iplot  = 1

#Nforced_batch_length = None ; # enforce the length of a batch (each batch will have a maximum of `Nforced_batch_length` records)
Nforced_batch_length = 2 ; # enforce the length of a batch (each batch will have a maximum of `Nforced_batch_length` records)
Nb_min_cnsctv = 2        ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)

min_nb_buoys_in_batch = 200 ; # minimum number of buoys for considering a batch a batch!

Nb_min_buoys = min_nb_buoys_in_batch ; # minimum number of buoys necessary to keep a given record of a given batch, when saving files and figures

l_drop_tooclose = False ; # PR: keep the one with the longest record...
NbPass = 2  # number of passes...
#
l_drop_overlap = True
#
rd_tol_km = 5.5 # tolerance distance in km below which we decide to cancel one of the 2 buoys! => for both `l_drop_tooclose` & `l_drop_overlap`

FillValue = -9999.

#================================================================================================


def __summary__( pNBini, pTcini, pIDs, pNRc ):
    (Nbtch,NbMax) = np.shape(pIDs)
    if Nbtch != len(pNBini):
        print('ERROR: [__summary__()] => error #1')
    if not NbMax == np.max(pNBini):
        print('ERROR: [__summary__()] => error #2')
    print('\n ==========   SUMMARY   ==========')
    print(' *** Number of identified batches: '+str(Nbtch))
    print(' *** Number of buoys selected in each batch:')
    for js in range(Nbtch):
        cTc0  = mjt.epoch2clock(pTcini[js])
        print('        * Batch #'+str(js)+' initiated at time bin centered around '+cTc0+' => has '+str(pNBini[js])+' buoys')
    print(' *** Max number of buoys possibly found in a batch = ',NbMax)
    print('     * shape of ZIDs =', np.shape(pIDs))
    print('     * shape of ZNRc =', np.shape(pNRc))
    print(' ===================================\n')


#================================================================================================

if __name__ == '__main__':

    for cd in ['npz','nc','figs']:
        if not path.exists('./'+cd): mkdir('./'+cd)
    if not path.exists('./figs/SELECTION'): mkdir('./figs/SELECTION')

    ####################################################################################################
    narg = len(argv)
    if not narg in [5]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YYYYMMDD(_HH:MM)1> <YYYYMMDD(_HH:MM)2> <dt_binning (hours)>')
        exit(0)
    cf_in    =     argv[1]
    cdate1   =     argv[2]
    cdate2   =     argv[3]
    idtbin_h = int(argv[4])
    ####################################################################################################
    dt_bin_sec =   float(idtbin_h*3600) ; # bin width for time scanning in [s], aka time increment while
    #                                     # scanning for valid etime intervals

    cdt1, cdt2, cdtS1, cdtS2 = mjt.DateString( cdate1, cdate2, returnShort=True )
    
    # File to read selected batches (generated by `00_scan_batches.py`) from:
    cf_npz_in = './npz/RGPS_batch_selection_'+cdtS1+'_'+cdtS2+'.npz'
    mjt.chck4f(cf_npz_in)

    max_t_dev_allowed_in_bin = dt_bin_sec/2.01 ; # Inside a given time bin of a given batch, a point should not be further in time
    #                                           # to the time mean of all points of this time bin than `max_t_dev_allowed_in_bin`

    if l_drop_overlap and l_drop_tooclose:
        print(' ERROR: you cannot use `l_drop_overlap` and `l_drop_tooclose`! Choose one of the two!!!')
        exit(0)

    print('\n *** Date range to restrain data to:')
    print(' ==> '+cdt1+' to '+cdt2 )

    idt1, idt2 = mjt.clock2epoch(cdt1), mjt.clock2epoch(cdt2)
    print( '   ===> in epoch time: ', idt1, 'to', idt2 )
    print( '       ====> double check: ', mjt.epoch2clock(idt1), 'to',  mjt.epoch2clock(idt2))

    if Nforced_batch_length:
        if Nb_min_cnsctv > Nforced_batch_length:
            print('ERROR: `Nb_min_cnsctv` cannot be > `Nforced_batch_length` !'); exit(0)

    # Build scan time axis willingly at relative high frequency (dt_bin_sec << dt_buoy_Nmnl)
    NTbin, vTbin = mjt.TimeBins4Scanning( idt1, idt2, dt_bin_sec, iverbose=idebug-1 )


    # Load data prepared for the time-range of interest (arrays are masked outside, except vtime0!)
    Np, Nb, vIDsWP, vtime0, vIDs0, vlat0, vlon0, vykm0, vxkm0 = mjt.LoadData4TimeRange( idt1, idt2, cf_in, l_doYX=True )
    # * Np: number of points of interst
    # * Nb: number of unique buoys of interest
    # * vIDsWP: unique IDs of the buoys of interest



    print('\n *** Opening "selected batches" data into '+cf_npz_in+'!')

    with np.load(cf_npz_in) as data:
        Nbatches = data['Nbatches']
        ZNB_ini  = data['ZNB_ini']
        ZTc_ini  = data['ZTc_ini']
        ZIDs     = data['IDs']
        ZNRc     = data['NRc']
        ZIX0     = data['ZIX0']
    # For some reason, masked shit not preserved via savez/load...
    ZIDs = np.ma.masked_where( ZIDs<0, ZIDs )
    ZNRc = np.ma.masked_where( ZNRc<0, ZNRc )
    #
    print('      => we have '+str(Nbatches)+' batches!')

    
    Ncsrec_max = np.max(ZNRc) ; # maximum number of valid consecutive records for a buoy
    for jr in range(Ncsrec_max):
        ZIX0[:,:,jr] = np.ma.masked_where( ZIX0[:,:,jr]==-999, ZIX0[:,:,jr] ) ; # fixme: check!


    # Reminder of what we found with previous script:
    mjt.batchSummaryRGPS( ZNB_ini, ZTc_ini, ZIDs, ZNRc )


    # Loop along batches:
    for jS in range(Nbatches):

        cs   = str(jS)
        vIDs = np.ma.MaskedArray.compressed( ZIDs[jS,:] ) ; # valid IDs for current batch: shrinked, getting rid of masked points
        NvB  = ZNB_ini[jS]
        rTc  = ZTc_ini[jS]
        cTc  = mjt.epoch2clock(rTc)

        if NvB != len(vIDs): print('ERROR Z1!'); exit(0)
        print('\n *** Having a look at batch #'+cs+' initiated for time bin centered around '+cTc+' !')

        if Nforced_batch_length:
            NCRmax = Nforced_batch_length
        else:
            NCRmax = np.max(ZNRc[jS,:]) ; # Max number of record from the buoy that has the most
            print('     ===> has '+str(NvB)+' valid buoys at start!')
            print('     ===> the buoy with most records has '+str(NCRmax)+' of them!')

        if idebug>2:
            print('        => with following IDs:')
            for jii in vIDs: print(jii,' ', end='')
            print('')


        # Creating arrays specific to a single batch to save (and plot from)
        #####################################################################

        nBpR = np.zeros( NCRmax      , dtype=int) ; # number of remaining buoys at given record
        xXkm = np.zeros((NCRmax,NvB)) + FillValue
        xYkm = np.zeros((NCRmax,NvB)) + FillValue
        xlon = np.zeros((NCRmax,NvB)) + FillValue
        xlat = np.zeros((NCRmax,NvB)) + FillValue
        xtim = np.zeros((NCRmax,NvB)) + FillValue      ; # the exact time for each buoy!
        xmsk = np.zeros((NCRmax,NvB) , dtype='i1')  ; # the mask for exluding buoys that stick out in time...
        xix0 = np.zeros((NCRmax,NvB) , dtype=int)  ;

        for jb in range(NvB):
            #
            nvr = ZNRc[jS,jb] ; # how many successive valid records for this buoy (at least `Nb_min_cnsctv`)
            if Nforced_batch_length:
                nvr = min( nvr, Nforced_batch_length )
            #
            if nvr<Nb_min_cnsctv: print('ERROR Z2!'); exit(0)
            #
            xix0[0:nvr,jb] = ZIX0[jS,jb,0:nvr] ; #lolo # all consecutive point position (indices as in `*0` arrays) for thi buoy
            #
            #(idx0_id,) = np.where( vIDs0 == vIDs[jb])
            #indv = idx0_id[0:nvr] ; # from record `nvr` onward buoy has been canceled (due to rogue time / expected time)
            indv = xix0[0:nvr,jb].copy() #
            #
            xXkm[0:nvr,jb] = vxkm0[indv]
            xYkm[0:nvr,jb] = vykm0[indv]
            xmsk[0:nvr,jb] = 1
            xlon[0:nvr,jb] = vlon0[indv]
            xlat[0:nvr,jb] = vlat0[indv]
            xtim[0:nvr,jb] = vtime0[indv]


        nBpR[:] = [ np.sum(xmsk[jr,:]) for jr in range(NCRmax) ] ; # How many buoys still present at each record?
        if np.max(nBpR) != NvB: print('ERROR: max(nBpR) != NvB !'); exit(0)
        if idebug>1: print('     +++ num. of boys still present at each record of batch #'+cs+':',nBpR[:])

        # A "mean" time axis along records based on mean accross the buoys remaining at this record...
        ztim = xtim.copy()
        ztim = np.ma.masked_where( xmsk==0, ztim ) ; # otherwize the `mean` in next line would use zeros!!!!
        vtim = np.mean(ztim, axis=1) ; # average on the buoy axis, so `vtim` only dimension is records...

        # Nearest interpolation of vtim on the VTbin calendar !!!
        #    => so VT contains the mean date for all buoys at a given record but
        #       corresponding to a value taken from VTbin
        VT = np.zeros( (NCRmax,3), dtype=int )
        ir=0
        for rt in vtim:
            kp = np.argmin( np.abs(vTbin[:,0]-rt) )
            lBatchOk = ( rt>=vTbin[kp,1] and rt<vTbin[kp,2] )
            if not lBatchOk:
                print(' WARNING: could not locate mean buoy time `rt` in the identified time bin :() ')
                print('          => rt, and time bounds =', mjt.epoch2clock(rt), mjt.epoch2clock(vTbin[kp,1]), mjt.epoch2clock(vTbin[kp,2]))
                print('          => identified time bin is #'+str(kp+1)+' out of '+str(NTbin)+' in total')
                print('          ==> forget this batch!\n')
                break
            #
            if idebug>1: print('      rt =',rt,' => ',mjt.epoch2clock(rt),' => nearest of VTbin =',mjt.epoch2clock(vTbin[kp,0]))
            VT[ir,:] = vTbin[kp,:]
            ir+=1

        if lBatchOk:
            # Now, in each record of the batch we should exclude buoys which time position is not inside the expected time bin
            # or is just too far away from the mean of all buoys
            # => if such a buoy is canceld at batch # k, it should also be canceled at following records
            iFU, xmsk, nBpR = mjt.BatchTimeSanityCheck( cs, ztim, VT, xmsk, nBpR, max_t_dev_allowed_in_bin, iverbose=2 )
            #
            del ztim
    
            if iFU>0:
                #print('old shape =', np.shape(xmsk))
                # => we masked some first and/or second buoy records, so we can shrink the arrays accordingly
                (idxK0,) , (idxK1,) = np.where(xmsk[0,:]==1) , np.where(xmsk[1,:]==1)
                
                if len(idxK0)==0 or len(idxK1)==0:
                    print('FIXME: `len(idxK0)==0 or len(idxK1)==0`')
                    exit(0)
                
                idxK = np.unique( np.concatenate([idxK0,idxK1]) )
                NvB  = len(idxK)
                xmsk = xmsk[:,idxK]
                xXkm = xXkm[:,idxK]
                xYkm = xYkm[:,idxK]
                xlon = xlon[:,idxK]
                xlat = xlat[:,idxK]
                xtim = xtim[:,idxK]
                vIDs = vIDs[  idxK]
                #
                ztim = xtim.copy()
                ztim = np.ma.masked_where( xmsk==0, ztim ) ; # otherwize the `mean` in next line would use zeros!!!!
                vtim = np.mean(ztim, axis=1) ; # average on the buoy axis, so `vtim` only dimension is records...
                #
                #print('new shape =', np.shape(xmsk))
                
    
            if l_drop_tooclose:
                jr = 0 ; # we work with first record !!!
                #
                NvB, idxK = mjt.CancelTooClose( jr, rd_tol_km, xlat, xlon, xmsk, NbPass=2 )
                #
                xmsk = xmsk[jr:,idxK]
                xlat = xlat[jr:,idxK]
                xlon = xlon[jr:,idxK]
                xYkm = xYkm[jr:,idxK]
                xXkm = xXkm[jr:,idxK]
                xtim = xtim[jr:,idxK]
                vIDs =     vIDs[idxK]
                del idxK
                print('\n *** UPDATE: based on "almost-overlap" cleaning at scale of '+str(rd_tol_km)+' km => '+str(NvB)+' buoys left to follow!')
    
    
            # The one to keep!!!
            if l_drop_overlap:
                jr = 0 ; # we work with first record !!!
                (idxn,) = np.where( xmsk[jr,:]==1 )
                zcoor = np.array([ xXkm[jr,:], xYkm[jr,:] ]).T ; # for `gudhi` coord = [X,Y] !
                _, _, idx_keep = mjt.SubSampCloud( rd_tol_km, zcoor )
                idx_rm = np.setdiff1d( idxn, idx_keep ) ; # keep values of `idxn` that are not in `idx_keep`
                xmsk[jr:,idx_rm] = 0 ; # supressing at this time records and all those following!!!
                Nrm = len(idx_rm)
                #
                (idxK,) = np.where(xmsk[0,:]==1)
                NvB  = len(idxK)
                xmsk = xmsk[jr:,idxK]
                xlat = xlat[jr:,idxK]
                xlon = xlon[jr:,idxK]
                xYkm = xYkm[jr:,idxK]
                xXkm = xXkm[jr:,idxK]
                xtim = xtim[jr:,idxK]
                vIDs =     vIDs[idxK]
                del zcoor, idxn, idx_keep, idx_rm, idxK
                #
                print('\n *** UPDATE: "sub-sampling" cleaning at scale of '+str(rd_tol_km)+' km => '+str(NvB)+' buoys left to follow! ('+str(Nrm)+' removed)')
    
    
    
            ##############################
            # Time to save the stuff !!! #
            ##############################
    
    
            # Masking:
            xXkm = np.ma.masked_where( xmsk==0, xXkm )
            xYkm = np.ma.masked_where( xmsk==0, xYkm )
            xlon = np.ma.masked_where( xmsk==0, xlon )
            xlat = np.ma.masked_where( xmsk==0, xlat )
            xtim = np.ma.masked_where( xmsk==0, xtim )
            vIDs = np.ma.masked_where( xmsk[0,:]==0, vIDs )
    
            cout_root = 'SELECTION_RGPS_S'+'%3.3i'%(jS)
    
            if iplot>0:
                # Batch time evolution on Arctic map:
                kf = mjt.ShowBuoysMap_Trec( vtim, xlon, xlat, pvIDs=vIDs, cnmfig='SELECTION/'+cout_root,
                                            clock_res='d', NminPnts=Nb_min_buoys )
    
            # GENERATION OF COMPREHENSIVE NETCDF FILE:
            #  * 1 file per batch
            #  * for the time variable inside netCDF, we chose the mean time accross buoys in the bin used aka `vtim`
            #  * for the file name we chose the center of the bin used aka `VT[:,0]`
            cdt1, cdt2 = split(':',mjt.epoch2clock(VT[ 0,0]))[0] , split(':',mjt.epoch2clock(VT[-1,0]))[0]  ; # keeps at the hour precision...
            cdt1, cdt2 = str.replace( cdt1, '-', '') , str.replace( cdt2, '-', '')
            cdt1, cdt2 = str.replace( cdt1, '_', 'h') , str.replace( cdt2, '_', 'h')
            cf_nc_out = './nc/'+cout_root+'_'+cdt1+'_'+cdt2+'.nc'
    
            print('   * Batch =',jS,' => saving '+cf_nc_out)
            kk = mjt.ncSaveCloudBuoys( cf_nc_out, vtim, vIDs, xYkm, xXkm, xlat, xlon, mask=xmsk,
                                       xtime=xtim, fillVal=FillValue, corigin='RGPS' )
            

    ### for jS in range(Nbatches)
