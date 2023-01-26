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
#
# TO DO:
#  finish the inclusion of `XIX0` !!!
#
# Date dans les noms de fichier e figures geo completement naze!!! En dessous de l\intitialisation du stream!!!
#  => doit y a voir une sacree merde!
#  => la date de l'intitialisation du stream et celle de ses records suivant doit utiliser soit:
#     => le centre du bin (marche que pour 1er record)
#     => la date moyenne prise sur les bouees utilisee dans le record particulier d'un stream donne!!!!
#
# Faire un effort sur la date moyenne (centree?) du stream, et pour les raisons suivante:
# - quel bouee doublon exclure si plusieurs memes ID de boue dans le meme bin
# - etre precis dans les noms de fichier et figures!!!
# - savoir quand on peut recomposer des streams avec un des scrips suivants

# * find the reason why the stream willing to start at `1997-01-04_18h00` will:
#   - work for short selection (1st to 10th of January)
#   - be cancelled for longer selection (1st to 31st of January)
# * also make the reason(s) why a stream was canceled clear, they could be many, will help!
#

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np

from re import split
from scipy import interpolate
from climporn import epoch2clock, clock2epoch
import mojito as mjt

idebug = 1

cdt_pattern = 'YYYY-MM-DD_00:00:00' ; # pattern for dates

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

dt_buoy_Nmnl = 3*24*3600     ; # the expected nominal time step of the input data, ~ 3 days [s]
max_dev_from_dt_buoy_Nmnl =  12.*3600. ; # maximum allowed deviation from the `dt_buoy_Nmnl` between 2 consecutive records of buoy [s]

Ns_max  =  100 ; # Max number of Streams, guess!, just for dimensionning array before knowing!!!
NrB_max =  50  ; # Max number of valid consecutive records for a given buoy, guess!, just for dimensionning array before knowing!!!

min_nb_buoys_in_stream = 200 ; # minimum number of buoys for considering a stream a stream!

Nb_min_cnsctv = 2   ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)

MinDistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]

Nb_min_buoys = min_nb_buoys_in_stream ; # minimum number of buoys necessary to keep a given record of a given stream, when saving files and figures

list_expected_var = [ 'index', 'x', 'y', 'lon', 'lat', 'q_flag', 'time' ]

interp_1d = 0 ; # Time interpolation to fixed time axis: 0 => linear / 1 => akima




def __summary__( pNBini, pTcini, pIDs, pNRc ):
    (Nstrm,NbMax) = np.shape(pIDs)
    if Nstrm != len(pNBini):
        print('ERROR: [__summary__()] => error #1')
    if not NbMax == np.max(pNBini):
        print('ERROR: [__summary__()] => error #2')
    print('\n ==========   SUMMARY   ==========')
    print(' *** Number of identified streams: '+str(Nstrm))
    print(' *** Number of buoys selected in each stream:')
    for js in range(Nstrm):
        cTc0  = epoch2clock(pTcini[js])
        print('        * Stream #'+str(js)+' initiated at time bin centered around '+cTc0+' => has '+str(pNBini[js])+' buoys')
    print(' *** Max number of buoys possibly found in a stream = ',NbMax)
    print('     * shape of ZIDs =', np.shape(pIDs))
    print('     * shape of ZNRc =', np.shape(pNRc))
    print(' ===================================\n')



if __name__ == '__main__':

    cdata_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
    fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    for cd in ["npz","figs"]:
        if not path.exists('./'+cd): mkdir('./'+cd)
    if not path.exists('./figs/SELECTION'): mkdir('./figs/SELECTION')

    ####################################################################################################
    narg = len(argv)
    if not narg in [5]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YYYYMMDD1> <YYYYMMDD2> <dt_binning (hours)>')
        exit(0)
    cf_in    =     argv[1]
    cYmmd1   =     argv[2]
    cYmmd2   =     argv[3]
    idtbin_h = int(argv[4])
    ####################################################################################################
    
    dt_bin =   float(idtbin_h*3600) ; # bin width for time scanning in [s], aka time increment while scanning for valid time intervals

    cY1,  cY2  = cYmmd1[0:4], cYmmd2[0:4]
    cmm1, cmm2 = cYmmd1[4:6], cYmmd2[4:6]
    cdd1, cdd2 = cYmmd1[6:8], cYmmd2[6:8]

    cdt1 = str.replace(cdt_pattern,'YYYY',cY1) ; cdt1 = str.replace(cdt1,'MM',cmm1) ; cdt1 = str.replace(cdt1,'DD',cdd1)
    cdt2 = str.replace(cdt_pattern,'YYYY',cY2) ; cdt2 = str.replace(cdt2,'MM',cmm2) ; cdt2 = str.replace(cdt2,'DD',cdd2)

    cf_out = 'RGPS_ice_drift_'+split('_', cdt1)[0]+'_'+split('_', cdt2)[0]+'_lb.nc' ;# netCDF file to generate

    # File to save work at intermediate stage
    cf_npz_intrmdt = './npz/'+str.replace( path.basename(cf_in), '.nc4', '.npz' )


    print("\n *** Date range to restrain data to:")
    print(" ==> "+cdt1+" to "+cdt2 )

    rdt1, rdt2 = clock2epoch(cdt1), clock2epoch(cdt2)
    print( "   ===> in epoch time: ", rdt1, "to", rdt2 )
    print( "       ====> double check: ", epoch2clock(rdt1), "to",  epoch2clock(rdt2))

    # Load `distance to coast` data:
    vlon_dist, vlat_dist, xdist = mjt.LoadDist2CoastNC( fdist2coast_nc )

    # Build scan time axis willingly at relative high frequency (dt_bin << dt_buoy_Nmnl)
    NTbin, vTbin, cTbin =   mjt.TimeBins4Scanning( rdt1, rdt2, dt_bin, iverbose=idebug-1 )
    
    # Open, inspect the input file and load raw data:
    Np0, vtime0, vx0, vy0, vlon0, vlat0, vBIDs0 = mjt.InspectLoadData( cf_in, list_expected_var )
    
    # Masking all point that are before and beyond our period of interest (note: `vtime0` is not masked!):
    Nb, vIDsWP = mjt.KeepDataInterest( rdt1, rdt2, vtime0, vBIDs0, vx0, vy0, vlon0, vlat0,  rmskVal=-99999. )
    # * Nb: number of different buoys that exist for at least 1 record during specified date range aka whole period (WP)
    # * vIDsWP : array(Nb) list (unique) of IDs for these buoys
    
    if not path.exists(cf_npz_intrmdt):

        # Arrays along streams and buoys:
        # In the following, both Ns_max & Nb are excessive upper bound values....
        VTc_ini = np.zeros( Ns_max                ) - 999.; # time of the center of time bin that first detected this stream =~ initialization date for a given stream
        VNB_ini = np.zeros( Ns_max,      dtype=int) - 999 ; # number of valid buoys at initialization of each stream (normally max, as the number of buoys alive decreases as records go on...)
        XIDs    = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # stores buoys IDs in use in a given stream
        XNRc    = np.zeros((Ns_max, Nb), dtype=int) - 999 ; # stores the number of records for each buoy in a given stream
        Xmsk    = np.zeros((Ns_max, Nb), dtype=int)       ; # tells if given buoy of given stream is alive (1) or dead (0)
        XIX0   = np.zeros((Ns_max, Nb, NrB_max), dtype=int) - 999 ;

        IDX_in_use_G = []  ; # keeps memory of points (indices) that have already been used by previous streams

        istream        = -1
        for jt in range(NTbin):
            #
            rTc = vTbin[jt,0] ; # center of the current time bin
            rTa = vTbin[jt,1] ; # begining of time bin
            #
            print("\n *** Selecting point indices that exist at "+cTbin[jt]+" +-"+str(int(dt_bin/2./3600))+"h!")
            (idxOK0,) = np.where( np.abs( vtime0[:] - rTc ) < dt_bin/2.-0.1 ) ; # yes, removing 0.1 second to `dt_bin/2.`
            Nok0 = len(idxOK0)
            #
            if Nok0>0:
                #
                print('     => we have ',Nok0,'such indices!')                
                zIDsOK0 = vBIDs0[idxOK0]
                # If the width of the time bin is large enough, the same buoy can exist more than once in the same time bin:
                ziDsU, idxU = np.unique(zIDsOK0, return_index=True)
                idxOKU = idxOK0[idxU] ; # because `idxU` are indices in the `zIDsOK0` world, not in the original `vBIDs0` world...
                del idxU
                
                Nok1 = len(ziDsU)
                if Nok1 < Nok0:
                    # LOLO: dans les doublons qui suit il faut juste garder la boue la plus proche en temps du centre du bin !!!! #fixme
                    #       => pour l'instant celui des 2 (voir 3,4,..) qui reste (vi le `np.unique` du dessus doit être random)
                    # indices of the doublons:
                    idxD = np.setdiff1d( idxOK0, idxOKU, assume_unique=True ) ; # keep the values of `idxOK0` that are not in `idxOKU`
                    zIDsD = vBIDs0[idxD]
                    print('     => some buoys exist more than once in the current date range selection!')
                    print('        (keeping 1 unique occurence leads to a removal of ',Nok0-Nok1,' points!)')
                    print('       ==> these buoys are: ',zIDsD)
                #
                print('     => we have ',Nok1,' unique buoys in these ',Nok0,' selected indices!')
                del zIDsOK0, idxOK0

                # Exclude points if index has already been used:
                idxT  = np.setdiff1d( idxOKU, np.array(IDX_in_use_G)) ;#, assume_unique=True ) ; # keep values of `idxOKU` that are not in `IDX_in_use_G`
                vIDsT = vBIDs0[idxT] ; # the buoys IDs we work with
                Nok = len(vIDsT)
                #
                if idebug>0:
                    # Sanity check: if any of the buoys found here do not belong to the whole-period reference buoy list `vIDsWP`:
                    vOUT = np.setdiff1d( vIDsT, vIDsWP) ; # keep the values of `vIDsT` that are not in `vIDsWP`
                    if len(vOUT)!=0: print('ERROR: some buoy IDs involved in this date range bin are not refenced in `vIDsWP` !!!'); exit(0)
                    #
                    print('     => '+str(Nok)+' buoys are still in the game! ('+str(Nok1-Nok)+' buoys removed because their index is already in use...)')
                    
                Nbuoys_stream = 0
                IDX_in_use_l  = []  ; # keeps memory of buoys that are already been included, but only at this stream level
                
                if Nok >= min_nb_buoys_in_stream:
                    
                    istream   = istream+1 ; # that's a new stream
                    if idebug>0: print('    => this date range is potentially the first of stream #'+str(istream)+', with '+str(Nok)+' buoys!')

                    # Now, loop on all the points involved in this date range:
                    jb = -1              ; # buoy counter...
                    for jidx in idxT:
                        #
                        jID = vBIDs0[jidx]
                        #
                        # #fixme: I have the feeling that `IDX_in_use_l` is unecessary???
                        if (not jidx in IDX_in_use_G) and (not jidx in IDX_in_use_l):

                            jb = jb + 1

                            nbRecOK, idx0_id, vt1b = mjt.ValidCnsctvRecordsBuoy( rTa, jidx, jID, vtime0, vBIDs0, np.array(IDX_in_use_G),
                                                                                 dt_buoy_Nmnl, max_dev_from_dt_buoy_Nmnl )
                            # * nbRecOK : number of valid consecutive records for this buoy
                            # * idx0_id : array of location indices (in the raw data arrays) for these valid records of this buoy
                            # * vt1b    : array of dates associated with all these records [s]
                            #print('---lolo: after `mjt.ValidCnsctvRecordsBuoy()` should have fixed the potential exclusion there...')
                            #lolo: I thinks that's where we cancel too many buoys because a time gap means a kill???
                            #exit(0)
                            
                            # We want at least `Nb_min_cnsctv` consecutive records for the buoy:
                            if nbRecOK >= Nb_min_cnsctv:
                                # We want the buoy to be located at least `MinDistFromLand` km off the coast                                
                                it1 = idx0_id[0]    ; # initial position for the buoy: #fixme: control all time records?
                                rd_ini = mjt.Dist2Coast( vlon0[it1], vlat0[it1], vlon_dist, vlat_dist, xdist )
                                if rd_ini > MinDistFromLand:
                                    #IDX_in_use_l.append(jidx)                  ; # point for first record of buoy `jID`!
                                    for ii in idx0_id: IDX_in_use_l.append(ii) ; # points for following records of `jID
                                    if not (jidx in IDX_in_use_l): print('ERROR: `not (jidx in IDX_in_use_l)`!'); exit(0)
                                    #
                                    Nbuoys_stream = Nbuoys_stream + 1   ; # this is another valid buoy for this stream
                                    Xmsk[istream,jb] = 1                ; # flag for valid point
                                    XIDs[istream,jb] = jID              ; # keeps memory of select buoy
                                    XNRc[istream,jb] = nbRecOK          ; # keeps memory of n. of valid consec. records
                                    XIX0[istream,jb,:nbRecOK] = idx0_id[:]
                                    
                                ### if rd_ini > MinDistFromLand
                            ### if nbRecOK >= Nb_min_cnsctv
                        ### if (not jidx in IDX_in_use_G) and (not jidx in IDX_in_use_l)
                    ### for jidx in idxT

                    if Nbuoys_stream >= min_nb_buoys_in_stream:
                        print('   +++ CONFIRMED VALID STREAM #'+str(istream)+' +++ => retained '+str(Nbuoys_stream)+' buoys!')
                        VNB_ini[istream] = Nbuoys_stream
                        VTc_ini[istream] = rTc
                        # Only now can we register the points indices we used into `IDX_in_use_G`:
                        for jidx in IDX_in_use_l: IDX_in_use_G.append(jidx)
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

        # MESSAGE:
        # *** Selection of buoys that exist at 1997-01-31_00:00:00 +-3h!
        # => 0 buoys satisfy this!
        # Traceback (most recent call last):
        # File "/home/laurent/DEV/rgps/./01_selection_xy.py", line 339, in <module>
        # Zmsk[js,:NvB] = Xmsk[js,indOK]
        # ValueError: could not broadcast input array from shape (5143,) into shape (5120,)


        # Masking arrays:
        ZIDs = np.ma.masked_where( Zmsk==0, ZIDs )
        ZNRc = np.ma.masked_where( Zmsk==0, ZNRc )
        for jr in range(Ncsrec_max):
            ZIX0[:,:,jr] = np.ma.masked_where( Zmsk==0, ZIX0[:,:,jr] )
        del Zmsk

        # Visualize buoy IDs in each stream:
        #for js in range(Nstreams):
        #    print('Stream #'+str(js)+' => ZIDs[js,:] =')
        #    for jb in range(Nbuoys_max): print( ZIDs[js,jb],' ',end='')
        #    print('')
        #exit(0)

        __summary__(ZNB_ini, ZTc_ini, ZIDs, ZNRc)

        print('\n *** Saving intermediate data into '+cf_npz_intrmdt+'!')
        np.savez_compressed( cf_npz_intrmdt, Nstreams=Nstreams, ZNB_ini=ZNB_ini, ZTc_ini=ZTc_ini, IDs=ZIDs, NRc=ZNRc, ZIX0=ZIX0 )



    else:

        print('\n *** Opening intermediate data into '+cf_npz_intrmdt+'!')

        with np.load(cf_npz_intrmdt) as data:
            Nstreams = data['Nstreams']
            ZNB_ini  = data['ZNB_ini']
            ZTc_ini  = data['ZTc_ini']
            ZIDs     = data['IDs']
            ZNRc     = data['NRc']
            ZIX0     = data['ZIX0']
        # For some reason, masked shit not preserved via savez/load...
        ZIDs = np.ma.masked_where( ZIDs==-999, ZIDs )
        ZNRc = np.ma.masked_where( ZNRc==-999, ZNRc )
        #
        Ncsrec_max = np.max(ZNRc) ; # maximum number of valid consecutive records for a buoy
        for jr in range(Ncsrec_max):
            ZIX0[:,:,jr] = np.ma.masked_where( ZIX0[:,:,jr]==-999, ZIX0[:,:,jr] ) ; # fixme: check!

        
        __summary__(ZNB_ini, ZTc_ini, ZIDs, ZNRc)

    ### if not path.exists(cf_npz_intrmdt)



    ############################################################################################################


    # Loop along streams:
    for js in range(Nstreams):
        
        cs   = str(js)
        vids = np.ma.MaskedArray.compressed( ZIDs[js,:] ) ; # valid IDs for current stream: shrinked, getting rid of masked points
        NvB  = ZNB_ini[js]
        rTc  = ZTc_ini[js]
        cTc  = epoch2clock(rTc)

        if NvB != len(vids): print('ERROR Z1!'); exit(0)
        print('\n *** Having a look at stream #'+cs+' initiated for time bin centered around '+cTc+' !')

        NCRmax = np.max(ZNRc[js,:]) ; # Max number of record from the buoy that has the most
        print('     ===> has '+str(NvB)+' valid buoys at start!')
        print('     ===> the buoy with most records has '+str(NCRmax)+' of them!')

        if idebug>2:
            print('        => with following IDs:')
            for jii in vids: print(jii,' ', end="")
            print('')


        # Creating arrays specific to a single stream to save (and plot from)
        #####################################################################

        nBpR = np.zeros( NCRmax      , dtype=int) ; # number of remaining buoys at given record
        xx   = np.zeros((NCRmax,NvB)) - 9999.
        xy   = np.zeros((NCRmax,NvB)) - 9999.
        xlon = np.zeros((NCRmax,NvB)) - 9999.
        xlat = np.zeros((NCRmax,NvB)) - 9999.
        xtim = np.zeros((NCRmax,NvB)) - 9999.      ; # the exact time for each buoy!
        xmsk = np.zeros((NCRmax,NvB) , dtype=int)  ; # the mask for exluding buoys that stick out in time...
        xix0 = np.zeros((NCRmax,NvB) , dtype=int)  ;

        for jb in range(NvB):
            #
            nvr = ZNRc[js,jb] ; # how many successive valid records for this buoy (at least `Nb_min_cnsctv`)
            if nvr<Nb_min_cnsctv: print('ERROR Z2!'); exit(0)
            #
            xix0[0:nvr,jb] = ZIX0[js,jb,0:nvr] ; #lolo # all consecutive point position (indices as in `*0` arrays) for thi buoy
            #
            #(idx0_id,) = np.where( vBIDs0 == vids[jb])
            #indv = idx0_id[0:nvr] ; # from record `nvr` onward buoy has been canceled (due to rogue time / expected time)
            indv = xix0[0:nvr,jb].copy()
            #
            xx[0:nvr,jb]   =    vx0[indv]
            xy[0:nvr,jb]   =    vy0[indv]
            xmsk[0:nvr,jb] =  1
            xlon[0:nvr,jb] =  vlon0[indv]
            xlat[0:nvr,jb] =  vlat0[indv]
            xtim[0:nvr,jb] = vtime0[indv]


        nBpR[:] = [ np.sum(xmsk[jr,:]) for jr in range(NCRmax) ] ; # How many buoys still present at each record?
        if np.max(nBpR) != NvB: print('ERROR: max(nBpR) != NvB !'); exit(0)
        if idebug>1: print('     +++ num. of boys still present at each record of stream #'+cs+':',nBpR[:])

        # There might be doublons in coordinates!
        ifd = 0
        for jr in range(NCRmax):
            xcoor = np.array([ xx[jr,:], xy[jr,:] ]).T
            idx_clones = mjt.idx_suppress_xy_copies( xcoor, rmask_val=-9999. )
            if len(idx_clones) > 0:
                # There are doublons!
                ifd = ifd + len(idx_clones) ; # another hit
                xmsk[jr,idx_clones] = 0
            del xcoor, idx_clones
        if ifd>0:
            # Again with the nlen(idx_clones)ew mask if we found doublon coordinates:
            nBpR[:] = [ np.sum(xmsk[jr,:]) for jr in range(NCRmax) ] ; # How many buoys still present at each record?
            NvB = np.max(nBpR)
            if idebug>1:
                print('     +++ we removed '+str(ifd)+' buoy records due to doublon of coordinates!')
                print('     +++ UPDATE: num. of boys still present at each record of stream #'+cs+':',nBpR[:])


        # A "mean" time axis for all the buoys
        ztim = xtim.copy()
        ztim = np.ma.masked_where( xmsk==0, ztim ) ; # otherwize the `mean` in next line would use zeros!!!!
        vtim = np.mean(ztim, axis=1) ; # average on the buoy axis, so `vtim` only dimension is records...
        del ztim

        # Nearest interpolation of vtim on the VTbin calendar !!!
        #    => so VT contains the mean date for all buoys at a given record but
        #       corresponding to a value taken from VTbin
        VT = vtim.copy()
        i=0
        for rt in vtim:
            idx = np.argmin( np.abs(vTbin[:,0]-rt) )
            if idebug>0: print('      rt =',rt,' => ',epoch2clock(rt),' => nearest of VTbin =',epoch2clock(vTbin[idx,0]))
            VT[i] = vTbin[idx,0]
            i=i+1

        # Are there buoys which are too far from this reference time ? #fixme
        #for jb in range(NvB):
        #    idt = np.abs(xtim[:,jb] - VT[:])
        #    if np.any(idt>dt_bin/2.):
        #        print('WOW, buoy #'+str(vids[jb])+' is more than '+str(int(dt_bin/2./3600))
        #              +' hours away from reference time...')
        #        #(idw,) = np.where(idt>dt_bin/2.)
        #        #print(idw)
        #        #print(' ==> supressing '+str(len(idw))+' values!')
        #        #xmsk[idw,jb] = 0


        if idebug>0:
            # Stream time evolution on Arctic map:
            kf = mjt.ShowBuoysMap_Trec( vtim, xlon, xlat, pvIDs=[], cnmfig='SELECTION/geo_buoys_RGPS_stream'+'%3.3i'%(js),
                                        clock_res='d', NminPnts=Nb_min_buoys )
        del vtim

        # Saving 1 file per stream and per record:
        for jr in range(NCRmax):

            Nbuoys = nBpR[jr] ; # number of buoys alive in current record of this stream
            idate =         int(VT[jr])
            cdate = epoch2clock(VT[jr])

            if Nbuoys >= Nb_min_buoys:
                ctr = str.replace( str.replace(cdate[0:16],':','h') ,'-','' ) ; # precision to the minute without ':'
                cf_out = './npz/SELECTION_buoys_RGPS_stream'+'%3.3i'%(js)+'_'+ctr+'.npz'

                vmsk = xmsk[jr,:]
                (indV,) = np.where(vmsk==1) ; # index of valid points
                if len(indV)!= Nbuoys: print('ERROR: len(indV)!= Nbuoys) !'); exit(0)

                np.savez_compressed( cf_out, itime=idate, date=cdate, Npoints=Nbuoys, vids=vids[indV],
                                     vtime=xtim[jr,indV], vx=xx[jr,indV], vy=xy[jr,indV], vlon=xlon[jr,indV], vlat=xlat[jr,indV] )

                if idebug>1:
                    # Plot on cartesian coordinates (km):
                    cfpng = './figs/SELECTION/xy_buoys_RGPS_stream'+'%3.3i'%(js)+'_'+ctr+'.png'
                    if jr==0:
                        zrx = [ np.min(xx[jr,indV])-100. , np.max(xx[jr,indV])+100. ]
                        zry = [ np.min(xy[jr,indV])-100. , np.max(xy[jr,indV])+100. ]
                    kg = mjt.ShowTQMesh( xx[jr,indV], xy[jr,indV], ppntIDs=vids[indV], cfig=cfpng, lGeoCoor=False,
                                         zoom=5, rangeX=zrx, rangeY=zry )
            else:
                if idebug>0:
                    print('     ===> NOT saving record #'+str(jr)+' of stream #'+cs+
                          ' (unsufficient n. of buoys alive:',Nbuoys,')')



    ### for js in range(Nstreams)
