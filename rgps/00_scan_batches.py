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

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

dt_Nmnl         = 3*24*3600 ; # the expected nominal time step of the input data, ~ 3 days [s]
#                             # => attention must be paid to `max_dev_dt_Nmnl` later in the code...

NbtchMax  =  200 ; # Max number of Batches, guess!, just for dimensionning array before knowing!!!

min_nb_buoys_in_batch = 10 ; # minimum number of buoys for considering a batch a batch!

MinDistFromLand  = 100. ; # how far from the nearest coast should our buoys be? [km]

FillValue = -9999.


#================================================================================================

if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    for cd in ['npz','figs']:
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

    # Adjustements that depend on the width of bins:
    if dt_bin_sec < dt_Nmnl:
        max_dev_dt_Nmnl = 6*3600 ; # => 6 hours! maximum allowed deviation from the `dt_Nmnl` between 2 consecutive records of buoy [s]
    else:
        max_dev_dt_Nmnl = dt_Nmnl/3 ; # => ~ 1 day ! maximum allowed deviation from the `dt_Nmnl` between 2 consecutive records of buoy [s]

    print('\n *** Max. allowed deviation in time from the `dt_Nmnl` between 2 consec. points to select =',max_dev_dt_Nmnl/3600,'hours')
    
    cdt1, cdt2, cdtS1, cdtS2 = mjt.DateString( cdate1, cdate2, returnShort=True )
    
    # File to save work in:
    cf_npz_out = './npz/RGPS_batch_selection_'+cdtS1+'_'+cdtS2+'.npz'

    max_t_dev_allowed_in_bin = dt_bin_sec/2.01 ; # Inside a given time bin of a given batch, a point should not be further in time
    #                                           # to the time mean of all points of this time bin than `max_t_dev_allowed_in_bin`

    print('\n *** Date range to restrain data to:')
    print(' ==> '+cdt1+' to '+cdt2 )

    idt1, idt2 = mjt.clock2epoch(cdt1), mjt.clock2epoch(cdt2)
    print( '   ===> in epoch time: ', idt1, 'to', idt2 )
    print( '       ====> double check: ', mjt.epoch2clock(idt1), 'to',  mjt.epoch2clock(idt2))


    # Load `distance to coast` data:
    vlon_dist, vlat_dist, xdist = mjt.LoadDist2CoastNC( fdist2coast_nc )

    # Build binning time axis:
    NTbin, vTbin = mjt.TimeBins4Scanning( idt1, idt2, dt_bin_sec, iverbose=idebug-1 )


    # Load data prepared for the time-range of interest (arrays are masked outside, except vtime0!)
    Np, Nb, vIDsU0, vtime0, vIDs0, vlat0, vlon0 = mjt.LoadData4TimeRange( idt1, idt2, cf_in, l_doYX=False )
    # * Np: number of points of interst
    # * Nb: number of unique buoys of interest
    # * vIDsU0: unique IDs of the buoys of interest len=Nb  #fixme: sure?
    NpT = len(vtime0) ; # Actual length of the *0 arrays..
    
    # Arrays along batches and buoys:
    # In the following, both NbtchMax & Nb are excessive upper bound values....
    VTc_ini = np.zeros( NbtchMax                ) - 999.; # time at center of time bin that first detected this batch
    VNB_ini = np.zeros( NbtchMax,         dtype=int) - 999 ; # n. of valid buoys at inititialization of each batch
    XIDs    = np.zeros((NbtchMax, Nb),    dtype=int) - 999 ; # stores buoys IDs in use in a given batch
    XNRc    = np.zeros((NbtchMax, Nb),    dtype=int) - 999 ; # stores the number of records for each buoy in a given batch
    Xmsk    = np.zeros((NbtchMax, Nb),    dtype='i1')      ; # tells if given buoy of given batch is alive (1) or dead (0)
    XIX0    = np.zeros((NbtchMax, Nb, 2), dtype=int) - 999 ;

    IDXtakenG = []  ; # keeps memory of points (indices) that have already been used by previous batches

    ibatch = -1
    for jt in range(NTbin):
        #
        rTc = vTbin[jt,0] ; # center of the current time bin
        rTa = vTbin[jt,1] ; # begining of the current time bin
        rTb = vTbin[jt,2] ; # end of the current time bin
        #
        print('\n *** Selecting point pos. that exist at '+mjt.epoch2clock(rTc)+' +-'+str(int(dt_bin_sec/2./3600))+
              'h => between',mjt.epoch2clock(rTa),'&',mjt.epoch2clock(rTb) )

        (idxOK0,) = np.where( (vtime0>=rTa) & (vtime0<rTb) & (vIDs0.data>=0) )

        zIDsOK0 = vIDs0[idxOK0]
        ztimOK0 = vtime0[idxOK0]
        Nok0 = len(idxOK0)        
        print('     => after "inside time bin" selection: '+str(Nok0)+' pos. involving '+str(len(np.unique(zIDsOK0)))+' different buoys!')

        if Nok0>0:
            # If the width of the time bin is large enough (normally>3days),
            # the same buoy (ID) can have more than 1 position during the time bin,
            # => need to keep only one position:

            Nok0, idxOK0 = mjt.EMO2( zIDsOK0, ztimOK0, vIDs0, vtime0, idxOK0, rTc, criterion='nearest',
                                     dtNom=dt_Nmnl, devdtNom=max_dev_dt_Nmnl, iverbose=1 )
            del zIDsOK0, ztimOK0
            print('     => after "multi-occurence" exclusions: '+str(Nok0)+' pos. involving '+str(len(np.unique(vIDs0[idxOK0])))+' different buoys!')

            # Exclude points if index has already been used or canceled:
            idxOK  = np.setdiff1d( idxOK0, np.array(IDXtakenG)) ; # keep values of `idxOK0` that are not in `IDXtakenG`
            Nok    = len(idxOK)
            zIDsOK = vIDs0[idxOK] ; # the buoys IDs we work with
            if len(np.unique(zIDsOK)) != Nok:
                print('ERROR: `unique(zIDsOK) != Nok` => `ExcludeMultiOccurences()` did not do its job :('); exit(0)
            print('     => after "already in use" exclusions: '+str(Nok)+' pos. involving '+str(len(np.unique(zIDsOK)))+' different buoys!')

            if idebug>0:
                # Sanity check: if any of the buoys found here do not belong to the whole-period reference buoy list `vIDsU0`:
                vOUT = np.setdiff1d( zIDsOK, vIDsU0) ; # keep the values of `zIDsOK` that are not in `vIDsU0`
                if len(vOUT)!=0:
                    print('ERROR: the IDs of '+str(len(vOUT))+' buoys involved in this date range bin are not refenced in `vIDsU0` !!!')
                    print(' ==>', vOUT)
                    exit(0)
                print('     => '+str(Nok)+' buoys still in the game! ('+str(Nok0-Nok)+' removed because index already in use...)')

                
            #---------------------------------------------------------------------------------------------------------------------
            NBinStr  = 0     ; # number of buoys in the batch
            IDXofStr = []  ; # keeps memory of buoys that are already been included, but only at the batch level

            iBcnl_CR = 0  ; # counter for buoys excluded because of consecutive records...
            
            if Nok >= min_nb_buoys_in_batch:

                ibatch += 1 ; # that's a new batch
                if idebug>0: print('    => this date range is potentially the first of batch #'+str(ibatch)+', with '+str(Nok)+' buoys!')

                # Now, loop on all the remaining point positions involved in this time bin:
                jb = -1              ; # buoy index

                for jidx in idxOK:

                    jb += 1
                    
                    jID = vIDs0[jidx]

                    if jidx in IDXtakenG:
                        print('WOW! `jidx in IDXtakenG` !!!'); exit(0)
                    #if not jidx in IDXtakenG:

                    nbRecOK, idx0_id, vt1b = mjt.ValidNextRecord( rTa, jidx, vtime0, vIDs0, np.array(IDXtakenG),
                                                                  dt_Nmnl, max_dev_dt_Nmnl )

                    if nbRecOK==0:
                        IDXtakenG.append(jidx) ; # cancel jidx, it's the position of a mono-record buoy
                        #fixme: we should cancel this buoy GLOBALLY when nbRecOK==0, it's a mono-record buoy in the whole period of interest

                    # * nbRecOK : number of valid consecutive records for this buoy
                    # * idx0_id : array of location indices (in the raw data arrays) for these valid records of this buoy
                    # * vt1b    : array of dates associated with all these records [s]
                    
                    if nbRecOK == 2:
                        # We want the buoy to be located at least `MinDistFromLand` km off the coast
                        it1 = idx0_id[0]    ; # initial position for the buoy
                        rd_ini = mjt.Dist2Coast( vlon0[it1], vlat0[it1], vlon_dist, vlat_dist, xdist )
                        if rd_ini > MinDistFromLand:
                            IDXofStr.append(idx0_id[0]) ; # store point not to be used again. 0 because the 2nd record can be re-used!
                            #
                            NBinStr += 1   ; # this is another valid buoy for this batch
                            Xmsk[ibatch,jb] = 1                ; # flag for valid point
                            XIDs[ibatch,jb] = jID              ; # keeps memory of select buoy
                            XNRc[ibatch,jb] = nbRecOK          ; # keeps memory of n. of "retained" consec. records
                            XIX0[ibatch,jb,:nbRecOK] = idx0_id[:nbRecOK] ; # indices for these valid records of this buoy
                        ### if rd_ini > MinDistFromLand
                    else:
                        iBcnl_CR += 1
                    #
                    ### if nbRecOK == 2
                ### for jidx in idxOK
                print('     => '+str(iBcnl_CR)+' buoys were canceled for not having a reasonable upcomming position in time!')

                if NBinStr >= min_nb_buoys_in_batch:
                    print('   +++ C O N F I R M E D   V A L I D   B A T C H   #'+str(ibatch)+' +++ => selected '+str(NBinStr)+' buoys!')
                    VNB_ini[ibatch] = NBinStr
                    VTc_ini[ibatch] = rTc
                    # Only now can we register the points indices we used into `IDXtakenG`:
                    IDXtakenG.extend(IDXofStr)
                else:
                    print('  * Well, this batch did not make it through the selection process... :(')
                    Xmsk[ibatch,:] = 0
                    XIDs[ibatch,:] = -999 ; #rm ! masked later with Xmsk?
                    XNRc[ibatch,:] = -999 ; #rm !
                    XIX0[ibatch,:,:] = -999 ; #rm !

                    ibatch = ibatch - 1 ; # REWIND!
                    if idebug>0: print('    => this was not a batch! So back to batch #'+str(ibatch)+' !!!')

            ### if Nok > min_nb_buoys_in_batch

        else:
            print(' ==> no points to be found inside this period !!!')
        ### if Nok0>0
        print('')

    ### for jt in range(NTbin)


    Nbatches   = ibatch+1
    Nbuoys_max = np.max(VNB_ini)
    Ncsrec_max = np.max(XNRc) ; # maximum number of valid consecutive records for a buoy

    if Ncsrec_max != 2:
        print('ERROR: `Ncsrec_max != 2` !'); exit(0)

    # Now that we know how many batches and what is the maximum possible number of buoys into a batch,
    # we can reduce the arrays:
    ZTc_ini = np.zeros( Nbatches                        ) - 999.
    ZNB_ini = np.zeros( Nbatches             , dtype=int) - 999
    ZIDs    = np.zeros((Nbatches, Nbuoys_max), dtype=int) - 999 ; # bad max size!! Stores the IDs used for a given batch...
    ZNRc    = np.zeros((Nbatches, Nbuoys_max), dtype=int) - 999 ; # bad max size!! Stores the number of records
    Zmsk    = np.zeros((Nbatches, Nbuoys_max), dtype='i1')
    ZIX0    = np.zeros((Nbatches, Nbuoys_max, Ncsrec_max), dtype=int) - 999

    for js in range(Nbatches):
        (indOK,) = np.where(Xmsk[js,:]==1)
        NvB=VNB_ini[js]
        if len(indOK) != NvB:
            print('ERROR: len(indOK) != NvB !!!'); exit(0)
        #
        ZTc_ini[js]     = VTc_ini[js]
        ZNB_ini[js]     = VNB_ini[js]
        Zmsk[js,:NvB]   = Xmsk[js,indOK] ; # #fixme: crash in big run with message below:
        ZIDs[js,:NvB]   = XIDs[js,indOK]
        ZNRc[js,:NvB]   = XNRc[js,indOK]
        ZIX0[js,:NvB,:] = XIX0[js,indOK,:Ncsrec_max] ; #fixme okay?

    del Xmsk, VTc_ini, VNB_ini, XIDs, XNRc, XIX0

    # Masking arrays:
    ZIDs = np.ma.masked_where( Zmsk==0, ZIDs )
    ZNRc = np.ma.masked_where( Zmsk==0, ZNRc )
    for jr in range(Ncsrec_max):
        ZIX0[:,:,jr] = np.ma.masked_where( Zmsk==0, ZIX0[:,:,jr] )
    del Zmsk

    mjt.batchSummaryRGPS(ZNB_ini, ZTc_ini, ZIDs, ZNRc)

    print('\n *** Saving info about batches into: '+cf_npz_out+'!')
    np.savez_compressed( cf_npz_out, Nbatches=Nbatches, ZNB_ini=ZNB_ini, ZTc_ini=ZTc_ini, IDs=ZIDs, NRc=ZNRc, ZIX0=ZIX0 )

