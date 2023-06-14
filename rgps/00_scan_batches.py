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
from mojito import config as cfg

idebug = 0

NbtchMax  =  200 ; # Max number of Batches, guess!, just for dimensionning array before knowing!!!

#================================================================================================

if __name__ == '__main__':

    

    for cd in ['npz','figs']:
        if not path.exists('./'+cd): mkdir('./'+cd)
    if not path.exists('./figs/SELECTION'): mkdir('./figs/SELECTION')

    ####################################################################################################
    narg = len(argv)
    if not narg in [6]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YYYYMMDD(_HH:MM)1> <YYYYMMDD(_HH:MM)2> <dt_binning (hours)> <mode (rgps,model,xlose)>')
        exit(0)
    cf_in    =     argv[1]
    cdate1   =     argv[2]
    cdate2   =     argv[3]
    idtbin_h = int(argv[4])
    quality_mode = argv[5]

    ik = cfg.controlModeName( path.basename(__file__), quality_mode )

    ####################################################################################################

    kk = cfg.initialize( mode=quality_mode )
 
    dt_bin_sec =   float(idtbin_h*3600) ; # bin width for time scanning in [s], aka time increment while
    #                                     # scanning for valid etime intervals
    
    print('\n *** Max. allowed deviation in time from the `dt_Nmnl` between 2 consec. points, `rc_dev_dt_Nmnl` =',int(cfg.rc_dev_dt_Nmnl/3600),'h')

    cdt1, cdt2, cdtS1, cdtS2 = mjt.DateString( cdate1, cdate2, returnShort=True )

    # File to save work in:
    cf_npz_out = './npz/RGPS_batch_selection_dt'+str(idtbin_h)+'h_'+cdtS1+'_'+cdtS2+'_mode-'+quality_mode+'.npz'
    if path.exists(cf_npz_out):
        print('\n *** File '+cf_npz_out+' is already here!!! I have nothing to do!')
        exit(0)

    print('\n *** Date range to restrain data to:')
    print(' ==> '+cdt1+' to '+cdt2 )

    idt1, idt2 = mjt.clock2epoch(cdt1), mjt.clock2epoch(cdt2)
    print( '   ===> in epoch time: ', idt1, 'to', idt2 )
    print( '       ====> double check: ', mjt.epoch2clock(idt1), 'to',  mjt.epoch2clock(idt2))


    # Load `distance to coast` data:
    vlon_dist, vlat_dist, xdist = mjt.LoadDist2CoastNC( cfg.fdist2coast_nc )

    # Build binning time axis:
    NTbin, vTbin = mjt.TimeBins4Scanning( idt1, idt2, dt_bin_sec, iverbose=idebug-1 )


    # Load data prepared for the time-range of interest (arrays are masked outside, except vtime0!)
    Np, Nb, vIDsU0, vtime0, vIDs0, vlat0, vlon0 = mjt.LoadData4TimeRange( idt1, idt2, cf_in, l_doYX=False )
    # * Np: number of points of interst
    # * Nb: number of unique buoys of interest
    # * vIDsU0: unique IDs of the buoys of interest

    NpT = len(vtime0) ; # Actual length of the *0 arrays..

    # Arrays along batches and buoys:
    # In the following, both NbtchMax & Nb are excessive upper bound values....
    VjtBinN = np.zeros( NbtchMax     ,    dtype=int) - 999 ; # index that accesses the relevant time bin for the batch
    VNB_ini = np.zeros( NbtchMax,         dtype=int) - 999 ; # n. of valid buoys at inititialization of each batch
    XIDs    = np.zeros((NbtchMax, Nb),    dtype=int) - 999 ; # stores buoys IDs in use in a given batch
    XNRc    = np.zeros((NbtchMax, Nb),    dtype=int) - 999 ; # stores the number of records for each buoy in a given batch
    Xmsk    = np.zeros((NbtchMax, Nb),    dtype='i1')      ; # tells if given buoy of given batch is alive (1) or dead (0)
    XIX0    = np.zeros((NbtchMax, Nb, 2), dtype=int) - 999 ;
    #IDXtakenG = []  ; # keeps memory of buoy positions (indices) that have already been used by previous batches

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
                                     devdtNom=cfg.rc_dev_dt_Nmnl, iverbose=1 )
            del zIDsOK0, ztimOK0
            print('     => after "multi-occurence" exclusions: '+str(Nok0)+' pos. involving '+str(len(np.unique(vIDs0[idxOK0])))+' different buoys!')

            # Exclude points if index has already been used or canceled:
            #lolo:idxOK  = np.setdiff1d( idxOK0, np.array(IDXtakenG)) ; # keep values of `idxOK0` that are not in `IDXtakenG`
            idxOK = idxOK0;#lolo            
            Nok    = len(idxOK)
            zIDsOK = vIDs0[idxOK] ; # the buoys IDs we work with
            if len(np.unique(zIDsOK)) != Nok:
                print('ERROR: `unique(zIDsOK) != Nok` => `ExcludeMultiOccurences()` did not do its job :('); exit(0)
            #print('     => after "already in use" exclusions: '+str(Nok)+' pos. involving '+str(len(np.unique(zIDsOK)))+' different buoys!')

            if idebug>0:
                # Sanity check: if any of the buoys found here do not belong to the whole-period reference buoy list `vIDsU0`:
                vOUT = np.setdiff1d( zIDsOK, vIDsU0) ; # keep the values of `zIDsOK` that are not in `vIDsU0`
                if len(vOUT)!=0:
                    print('ERROR: the IDs of '+str(len(vOUT))+' buoys involved in this date range bin are not refenced in `vIDsU0` !!!')
                    print(' ==>', vOUT)
                    exit(0)
                print('     => '+str(Nok)+' buoys still in the game! ('+str(Nok0-Nok)+' removed because index already in use...)')


            #---------------------------------------------------------------------------------------------------------------------
            NBinBtch  = 0  ; # number of buoys in the batch
            iBcnl_CR  = 0  ; # counter for buoys excluded because of consecutive records...
            iBcnl_DC  = 0  ; # counter for buoys excluded because of excessive proximity to coastline

            if Nok >= cfg.nc_min_buoys_in_batch:

                ibatch += 1 ; # that's a new batch
                if idebug>0: print('    => this date range is potentially the first of batch #'+str(ibatch)+', with '+str(Nok)+' buoys!')

                # Now, loop on all the remaining point positions involved in this time bin:
                jb = -1              ; # buoy index
                for jidx in idxOK:

                    jb += 1

                    jID = vIDs0[jidx]

                    #lolo:
                    #if jidx in IDXtakenG:
                    #    print('WOW! `jidx in IDXtakenG` !!!'); exit(0) ; #fixme
                    ##if not jidx in IDXtakenG:

                    nbRecOK, idx0VUCR = mjt.ValidUpComingRecord( rTa, jidx, vtime0, vIDs0, devdtNom=cfg.rc_dev_dt_Nmnl )

                    #lolo:
                    #if nbRecOK==0:
                    #    IDXtakenG.append(jidx) ; # cancel jidx, it's the position of a mono-record buoy
                    #    #fixme: we should cancel this buoy GLOBALLY when nbRecOK==0, it's a mono-record buoy in the whole period of interest

                    # * nbRecOK : number of valid consecutive records for this buoy
                    # * idx0VUCR : array of location indices (in the raw data arrays) for these valid records of this buoy, including the present
                    #              `jidx` index at first position !

                    if nbRecOK == 2:
                        #
                        lOkDC = True
                        if cfg.lc_cancelShore:
                            # We want the buoy to be located at least `cfg.nc_MinDistFromLand` km off the coast
                            it1 = idx0VUCR[0]    ; # initial position for the buoy
                            rd_ini = mjt.Dist2Coast( vlon0[it1], vlat0[it1], vlon_dist, vlat_dist, xdist )
                            lOkDC = ( rd_ini > cfg.nc_MinDistFromLand )
                            del it1, rd_ini
                        #
                        if lOkDC:
                            NBinBtch += 1   ; # this is another valid buoy for this batch
                            Xmsk[ibatch,jb] = 1                ; # flag for valid point
                            XIDs[ibatch,jb] = jID              ; # keeps memory of select buoy
                            XNRc[ibatch,jb] = nbRecOK          ; # keeps memory of n. of "retained" consec. records
                            XIX0[ibatch,jb,:nbRecOK] = idx0VUCR[:nbRecOK] ; # indices for these valid records of this buoy
                        else:
                            iBcnl_DC += 1
                        #
                    else:
                        iBcnl_CR += 1
                    #
                    ### if nbRecOK == 2
                ### for jidx in idxOK
                
                if Nok-NBinBtch != iBcnl_CR+iBcnl_DC:
                    print('ERROR: `ok-NBinBtch != iBcnl_CR+iBcnl_DC` !!!',Nok-NBinBtch, iBcnl_CR+iBcnl_DC)
                    exit(0)
                #
                print('  *** '+str(Nok-NBinBtch)+' buoys were canceled:')
                print('         => '+str(iBcnl_CR)+' for not having a reasonable upcomming position in time!')
                if cfg.lc_cancelShore:
                    print('         => '+str(iBcnl_DC)+' for being excessively close to land! (d < '
                          +str(int(cfg.nc_MinDistFromLand))+'km)')
                else:
                    print('  *** Distance to coast is not a cause for concern as `lc_cancelShore==False` !')


                if NBinBtch >= cfg.nc_min_buoys_in_batch:
                    print('   +++ C O N F I R M E D   V A L I D   B A T C H   #'+str(ibatch)+' +++ => selected '+str(NBinBtch)+' buoys!')
                    VNB_ini[ibatch] = NBinBtch
                    VjtBinN[ibatch] = jt
                    # Only now can we register the points indices we used into `IDXtakenG`:
                    #lolo:IDXtakenG.extend(IDXofStr)
                else:
                    print('  * Well, this batch did not make it through the selection process... :(')
                    Xmsk[ibatch,:] = 0
                    XIDs[ibatch,:] = -999 ; #rm ! masked later with Xmsk?
                    XNRc[ibatch,:] = -999 ; #rm !
                    XIX0[ibatch,:,:] = -999 ; #rm !

                    ibatch = ibatch - 1 ; # REWIND!
                    if idebug>0: print('    => this was not a batch! So back to batch #'+str(ibatch)+' !!!')

            ### if Nok > cfg.nc_min_buoys_in_batch

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
    ZjtBinN = np.zeros( Nbatches             , dtype=int) - 999
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
        ZjtBinN[js]     = VjtBinN[js]
        ZNB_ini[js]     = VNB_ini[js]
        Zmsk[js,:NvB]   = Xmsk[js,indOK] ; # #fixme: crash in big run with message below:
        ZIDs[js,:NvB]   = XIDs[js,indOK]
        ZNRc[js,:NvB]   = XNRc[js,indOK]
        ZIX0[js,:NvB,:] = XIX0[js,indOK,:Ncsrec_max] ; #fixme okay?

    del Xmsk, VNB_ini, VjtBinN, XIDs, XNRc, XIX0

    # Masking arrays:
    ZIDs = np.ma.masked_where( Zmsk==0, ZIDs )
    ZNRc = np.ma.masked_where( Zmsk==0, ZNRc )
    for jr in range(Ncsrec_max):
        ZIX0[:,:,jr] = np.ma.masked_where( Zmsk==0, ZIX0[:,:,jr] )
    del Zmsk

    mjt.batchSummaryRGPS( vTbin, ZNB_ini, ZjtBinN, ZIDs, ZNRc)

    print('\n *** Saving info about batches into: '+cf_npz_out+'!')
    np.savez_compressed( cf_npz_out, dt_bin_sec=dt_bin_sec, vTbin=vTbin, Nbatches=Nbatches, NB_ini=ZNB_ini, jtBinN=ZjtBinN, IDs=ZIDs, NRc=ZNRc, ZIX0=ZIX0 )


