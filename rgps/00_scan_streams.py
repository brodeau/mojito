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

### TO DO: lili
# ExcludeMulitOccurences() should be able to use a criterion that cancels the point that has no reasonable successor (based on time step) !

# dangerous stuff is np.where(vIDs0>=0) because masked => can have the "--" value => better np.where(vIDs0.data>=0)


from sys import argv, exit
from os import path, environ, mkdir
import numpy as np

from re import split
from climporn import epoch2clock, clock2epoch
import mojito as mjt

idebug = 0

cdt_pattern = 'YYYY-MM-DD_hh:mm:00' ; # pattern for dates

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

#Nforced_stream_length = None ; # enforce the length of a stream (each stream will have a maximum of `Nforced_stream_length` records)
Nforced_stream_length = 2 ; # enforce the length of a stream (each stream will have a maximum of `Nforced_stream_length` records)
Nb_min_cnsctv = 2        ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)

dt_buoy_Nmnl = 3*24*3600     ; # the expected nominal time step of the input data, ~ 3 days [s]
max_dev_from_dt_buoy_Nmnl = 6*3600 ; # maximum allowed deviation from the `dt_buoy_Nmnl` between 2 consecutive records of buoy [s]

Ns_max  =  200 ; # Max number of Streams, guess!, just for dimensionning array before knowing!!!
NrB_max =  50  ; # Max number of valid consecutive records for a given buoy, guess!, just for dimensionning array before knowing!!!

min_nb_buoys_in_stream = 10 ; # minimum number of buoys for considering a stream a stream!

MinDistFromLand  = 100. ; # how far from the nearest coast should our buoys be? [km]

list_expected_var = [ 'index', 'x', 'y', 'lon', 'lat', 'time', 'idstream', 'streams' ]

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

    ldp1, lhp1 = len(cdate1)==8, len(cdate1)==14
    ldp2, lhp2 = len(cdate2)==8, len(cdate2)==14

    cY1,  cY2  = cdate1[0:4], cdate2[0:4]
    cmm1, cmm2 = cdate1[4:6], cdate2[4:6]
    cdd1, cdd2 = cdate1[6:8], cdate2[6:8]

    chhmm1, chhmm2 = '00:00', '00:00'
    if lhp1: chhmm1 = cdate1[9:11]+':'+cdate1[12:14]
    if lhp2: chhmm2 = cdate2[9:11]+':'+cdate2[12:14]

    cdt1 = str.replace(cdt_pattern,'YYYY-MM-DD',cY1+'-'+cmm1+'-'+cdd1)
    cdt2 = str.replace(cdt_pattern,'YYYY-MM-DD',cY2+'-'+cmm2+'-'+cdd2)
    cdt1 = str.replace(cdt1,'hh:mm',chhmm1)
    cdt2 = str.replace(cdt2,'hh:mm',chhmm2)

    # File to save work in:
    cf_npz_out = './npz/RGPS_stream_selection_'+cY1+cmm1+cdd1+'_'+cY2+cmm2+cdd2+'.npz'

    max_t_dev_allowed_in_bin = dt_bin_sec/2.01 ; # Inside a given time bin of a given stream, a point should not be further in time
    #                                           # to the time mean of all points of this time bin than `max_t_dev_allowed_in_bin`

    print('\n *** Date range to restrain data to:')
    print(' ==> '+cdt1+' to '+cdt2 )

    idt1, idt2 = clock2epoch(cdt1), clock2epoch(cdt2)
    print( '   ===> in epoch time: ', idt1, 'to', idt2 )
    print( '       ====> double check: ', epoch2clock(idt1), 'to',  epoch2clock(idt2))


    if Nforced_stream_length:
        if Nb_min_cnsctv > Nforced_stream_length:
            print('ERROR: `Nb_min_cnsctv` cannot be > `Nforced_stream_length` !'); exit(0)

    # Load `distance to coast` data:
    vlon_dist, vlat_dist, xdist = mjt.LoadDist2CoastNC( fdist2coast_nc )

    # Build scan time axis willingly at relative high frequency (dt_bin_sec << dt_buoy_Nmnl)
    NTbin, vTbin = mjt.TimeBins4Scanning( idt1, idt2, dt_bin_sec, iverbose=idebug-1 )


    # Load data prepared for the time-range of interest (arrays are masked outside, except vtime0!)
    Np, Nb, vIDsU0, vtime0, vIDs0, vlat0, vlon0 = mjt.LoadData4TimeRange( idt1, idt2, cf_in, list_expected_var, l_doYX=False )
    # * Np: number of points of interst
    # * Nb: number of unique buoys of interest
    # * vIDsU0: unique IDs of the buoys of interest len=Nb  #fixme: sure?
    NpT = len(vtime0) ; # Real length of the *0 arrays..
    
    #(idxM,) = np.where(vIDs0==145)
    #print('145 =>',len(idxM))
    #exit(0)

    
    # Exlude buoys and associate points corresponding to "mono-record" buoys (mono-record during the whole period, not during a bin)
    # - if the buoy ID related to a point exists only once in the whole period (not the bin),
    #   there is no reason to keep it
    #  => must 1st find the problematic buoys, and then cancel the associated point indices
    if 1==0:
        print(' Before mono-record stuff, we have '+str(Np),' points, for '+str(Nb)+' unique buoys.')
        (idxOK0,) = np.where(vIDs0>=0) ; # since vIDs0 is masked (time-range)
        if len(idxOK0) != Np:
            print('ERROR: `len(idxOK0) != Np`'); exit(0)
        #
        zmsk = np.zeros( NpT, dtype='i1' )
        zmsk[np.where(vIDs0.data>=0)] = 1    
        idxRM = []
        for kID in vIDsU0:
            (kdx,) = np.where( vIDs0== kID )
            ntp = len(kdx)
            if ntp < 2:
                #print(' !!! buoy with ID '+str(kID)+' has only 1 record in the period of interest !!!')
                if ntp != 1:
                    print('ERROR: `ntp != 1`'); exit(0)
                idxRM.append( kdx )
                zmsk[kdx] = 0
                
        idxK  = np.setdiff1d( idxOK0, np.array(idxRM)) ; # keep values of `idxOK0` that are not in `idxRM`
        #
        Np    = len(idxK)
        vIDsU0 = np.unique(vIDs0[idxK])
        Nb = len(vIDsU0)
        
        # Updating masked *0 arrays based on the mask: 
        vIDs0 = np.ma.masked_where( zmsk==0, vIDs0 )
        vlon0 = np.ma.masked_where( zmsk==0, vlon0 )
        vlat0 = np.ma.masked_where( zmsk==0, vlat0 )
        del idxRM, idxK, zmsk, idxOK0
        print('     => after "mono-record" buoys exclusions: '+str(Np)+' pos. involving '+str(Nb)+' different buoys!')

    
    # Arrays along streams and buoys:
    # In the following, both Ns_max & Nb are excessive upper bound values....
    VTc_ini = np.zeros( Ns_max                ) - 999.; # time at center of time bin that first detected this stream
    VNB_ini = np.zeros( Ns_max,      dtype=int) - 999 ; # n. of valid buoys at inititialization of each stream
    XIDs    = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # stores buoys IDs in use in a given stream
    XNRc    = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # stores the number of records for each buoy in a given stream
    Xmsk    = np.zeros((Ns_max, Nb), dtype=int)       ; # tells if given buoy of given stream is alive (1) or dead (0)
    XIX0    = np.zeros((Ns_max, Nb, NrB_max), dtype=int) - 999 ;

    IDXtakenG = []  ; # keeps memory of points (indices) that have already been used by previous streams

    istream        = -1
    for jt in range(NTbin):
        #
        rTc = vTbin[jt,0] ; # center of the current time bin
        rTa = vTbin[jt,1] ; # begining of the current time bin
        rTb = vTbin[jt,2] ; # end of the current time bin
        #
        print('\n *** Selecting point pos. that exist at '+epoch2clock(rTc)+' +-'+str(int(dt_bin_sec/2./3600))+
              'h => between',epoch2clock(rTa),'&',epoch2clock(rTb) )

        (idxOK0,) = np.where( (vtime0>=rTa) & (vtime0<rTb) & (vIDs0.data>=0) )

        zIDsOK0 = vIDs0[idxOK0]
        ztimOK0 = vtime0[idxOK0]
        Nok0 = len(idxOK0)        
        print('     => after "inside time bin" selection: '+str(Nok0)+' pos. involving '+str(len(np.unique(zIDsOK0)))+' different buoys!')

        if Nok0>0:
            # If the width of the time bin is large enough (normally>3days),
            # the same buoy ID can exist more than once in the same time bin,
            # and so also in `zIDsOK0`!
            if dt_bin_sec < dt_buoy_Nmnl:
                # Time bins are narrower than the nominal time step of the RGPS data...
                #   => we keep buoy occurence closest to that of center of current time bin
                Nok0, idxOK0 = mjt.ExcludeMulitOccurences( zIDsOK0, ztimOK0, vIDs0, idxOK0, rTc, criterion='center', iverbose=idebug )
            else:
                # Time bins are wider than the nominal time step of the RGPS data...
                #   => we keep buoy occurence with the earliest occurence
                Nok0, idxOK0 = mjt.ExcludeMulitOccurences( zIDsOK0, ztimOK0, vIDs0, idxOK0, rTc, criterion='first', iverbose=idebug )
            #
            del zIDsOK0, ztimOK0
            print('     => after "multi-occurence" exclusions: '+str(Nok0)+' pos. involving '+str(len(np.unique(vIDs0[idxOK0])))+' different buoys!')
            
            # Exclude points if index has already been used:
            idxOK  = np.setdiff1d( idxOK0, np.array(IDXtakenG)) ; # keep values of `idxOK0` that are not in `IDXtakenG`
            Nok    = len(idxOK)
            zIDsOK = vIDs0[idxOK] ; # the buoys IDs we work with
            if len(np.unique(zIDsOK)) != Nok:
                print('ERROR: `unique(zIDsOK) != Nok` => `ExcludeMulitOccurences()` did not do its job :('); exit(0)
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
            NBinStr  = 0     ; # number of buoys in the stream
            IDXofStr = []  ; # keeps memory of buoys that are already been included, but only at the stream level

            iBcnl_CR = 0  ; # counter for buoys excluded because of consecutive records...
            
            if Nok >= min_nb_buoys_in_stream:

                istream += 1 ; # that's a new stream
                if idebug>0: print('    => this date range is potentially the first of stream #'+str(istream)+', with '+str(Nok)+' buoys!')

                # Now, loop on all the remaining point positions involved in this time bin:
                jb = -1              ; # buoy index
                for jidx in idxOK:
                    #
                    jID = vIDs0[jidx]
                    #
                    if not jidx in IDXtakenG:

                        jb += 1

                        if Nforced_stream_length==2:
                            nbRecOK, idx0_id, vt1b = mjt.ValidNextRecord( rTa, jID, jidx, vtime0, vIDs0, np.array(IDXtakenG),
                                                                 dt_buoy_Nmnl, max_dev_from_dt_buoy_Nmnl )

                            if idebug>1:
                                if nbRecOK==0:
                                    print(' LOLO: this buoy is a mono-record buoy !!!')
                                    (idxM,) = np.where(vIDs0==jID)
                                    if len(idxM)!=1: print('ERROR ZXM!'); exit(0)
                                elif nbRecOK==1:
                                    print(' LOLO: not mono-record buoy but did not find a reasonable sucessor for it')
                                    (idxM,) = np.where((vIDs0==jID) & (vtime0>=rTa))
                                    print(' Dates for point + successors are:')
                                    for zt in vtime0[idxM]: print(epoch2clock(zt))
                                elif nbRecOK==2:
                                    print(' LOLO: WE FOUND a reasonable sucessor for this buoy')
                                    (idxM,) = np.where((vIDs0==jID) & (vtime0>=rTa))
                                    print(' Dates for point + all possible successors are:')
                                    for zt in vtime0[idxM]: print(epoch2clock(zt))
                                    print(' ==> the time selected is',epoch2clock(vtime0[idx0_id[1]]) )
                                    print(' ====> vt1b =', epoch2clock(vt1b[0]), epoch2clock(vt1b[1]))

                            #if nbRecOK==0:
                            #fixme: we should cancel this buoy GLOBALLY when nbRecOK==0, it's a mono-record buoy in the whole period of interest

                        else:                        
                            nbRecOK, idx0_id, vt1b = mjt.ValidCnsctvRecordsBuoy( rTa, jidx, vtime0, vIDs0, np.array(IDXtakenG),
                                                                                 dt_buoy_Nmnl, max_dev_from_dt_buoy_Nmnl )


                        # * nbRecOK : number of valid consecutive records for this buoy
                        # * idx0_id : array of location indices (in the raw data arrays) for these valid records of this buoy
                        # * vt1b    : array of dates associated with all these records [s]

                        
                        # We want at least `Nb_min_cnsctv` consecutive records for the buoy:
                        if nbRecOK >= Nb_min_cnsctv:
                            # We want the buoy to be located at least `MinDistFromLand` km off the coast
                            it1 = idx0_id[0]    ; # initial position for the buoy
                            rd_ini = mjt.Dist2Coast( vlon0[it1], vlat0[it1], vlon_dist, vlat_dist, xdist )
                            if rd_ini > MinDistFromLand:
                                IDXofStr.extend(idx0_id[:nbRecOK-1]) ; # store point not to be used again. -1 because the last record can be re-used!
                                #
                                NBinStr += 1   ; # this is another valid buoy for this stream
                                Xmsk[istream,jb] = 1                ; # flag for valid point
                                XIDs[istream,jb] = jID              ; # keeps memory of select buoy
                                XNRc[istream,jb] = nbRecOK          ; # keeps memory of n. of "retained" consec. records
                                XIX0[istream,jb,:nbRecOK] = idx0_id[:nbRecOK] ; # indices for these valid records of this buoy
                            ### if rd_ini > MinDistFromLand
                        else:
                            iBcnl_CR += 1
                        #
                        ### if nbRecOK >= Nb_min_cnsctv
                    ### if not jidx in IDXtakenG
                ### for jidx in idxOK
                print('     => '+str(iBcnl_CR)+' buoys were canceled for not having "proper" upcomming positions!')

                if NBinStr >= min_nb_buoys_in_stream:
                    print('   +++ C O N F I R M E D   V A L I D   S T R E A M   #'+str(istream)+' +++ => selected '+str(NBinStr)+' buoys!')
                    VNB_ini[istream] = NBinStr
                    VTc_ini[istream] = rTc
                    # Only now can we register the points indices we used into `IDXtakenG`:
                    IDXtakenG.extend(IDXofStr)
                else:
                    print('  * Well, this stream did not make it through the selection process... :(')
                    Xmsk[istream,:] = 0
                    XIDs[istream,:] = -999 ; #rm ! masked later with Xmsk?
                    XNRc[istream,:] = -999 ; #rm !
                    XIX0[istream,:,:] = -999 ; #rm !

                    istream = istream - 1 ; # REWIND!
                    if idebug>0: print('    => this was not a stream! So back to stream #'+str(istream)+' !!!')

            ### if Nok > min_nb_buoys_in_stream

        else:
            print(' ==> no points to be found inside this period !!!')
        ### if Nok0>0
        print('')

    ### for jt in range(NTbin)


    Nstreams   = istream+1
    Nbuoys_max = np.max(VNB_ini)
    Ncsrec_max = np.max(XNRc) ; # maximum number of valid consecutive records for a buoy


    # Now that we know how many streams and what is the maximum possible number of buoys into a stream,
    # we can reduce the arrays:
    ZTc_ini = np.zeros( Nstreams                        ) - 999.
    ZNB_ini = np.zeros( Nstreams             , dtype=int) - 999
    ZIDs    = np.zeros((Nstreams, Nbuoys_max), dtype=int) - 999 ; # bad max size!! Stores the IDs used for a given stream...
    ZNRc    = np.zeros((Nstreams, Nbuoys_max), dtype=int) - 999 ; # bad max size!! Stores the number of records
    Zmsk    = np.zeros((Nstreams, Nbuoys_max), dtype=int)
    ZIX0    = np.zeros((Nstreams, Nbuoys_max, Ncsrec_max), dtype=int) - 999

    for js in range(Nstreams):
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

    mjt.streamSummaryRGPS(ZNB_ini, ZTc_ini, ZIDs, ZNRc)

    print('\n *** Saving info about streams into: '+cf_npz_out+'!')
    np.savez_compressed( cf_npz_out, Nstreams=Nstreams, ZNB_ini=ZNB_ini, ZTc_ini=ZTc_ini, IDs=ZIDs, NRc=ZNRc, ZIX0=ZIX0 )


