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

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np

from re import split
from climporn import epoch2clock, clock2epoch
import mojito as mjt

idebug = 0
iplot  = 1

cdt_pattern = 'YYYY-MM-DD_hh:mm:00' ; # pattern for dates

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

#Nforced_stream_length = None ; # enforce the length of a stream (each stream will have a maximum of `Nforced_stream_length` records)
Nforced_stream_length = 2 ; # enforce the length of a stream (each stream will have a maximum of `Nforced_stream_length` records)
Nb_min_cnsctv = 2        ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)

min_nb_buoys_in_stream = 200 ; # minimum number of buoys for considering a stream a stream!

Nb_min_buoys = min_nb_buoys_in_stream ; # minimum number of buoys necessary to keep a given record of a given stream, when saving files and figures

list_expected_var = [ 'index', 'x', 'y', 'lon', 'lat', 'q_flag', 'time' ]

#================================================================================================

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
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YYYYMMDD(_HH:MM)1> <YYYYMMDD(_HH:MM)2> <dt_binning (hours)>')
        exit(0)
    cf_in    =     argv[1]
    cdate1   =     argv[2]
    cdate2   =     argv[3]
    idtbin_h = int(argv[4])
    ####################################################################################################        
    dt_bin_sec =   float(idtbin_h*3600) ; # bin width for time scanning in [s], aka time increment while
    #                                 # scanning for valid etime intervals
    #
    ldp1, lhp1 = len(cdate1)==8, len(cdate1)==14
    ldp2, lhp2 = len(cdate2)==8, len(cdate2)==14
    #
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
    
    cf_out = 'RGPS_ice_drift_'+split('_', cdt1)[0]+'_'+split('_', cdt2)[0]+'_lb.nc' ;# netCDF file to generate

    # File to save work at intermediate stage
    cf_npz_strinfo = './npz/'+str.replace( path.basename(cf_in), '.nc4', '.npz' )


    max_t_dev_allowed_in_bin = dt_bin_sec/2.01 ; # Inside a given time bin of a given stream, a point should not be further in time
    #                                           # to the time mean of all points of this time bin than `max_t_dev_allowed_in_bin`

    
    print("\n *** Date range to restrain data to:")
    print(" ==> "+cdt1+" to "+cdt2 )
    
    rdt1, rdt2 = clock2epoch(cdt1), clock2epoch(cdt2)
    print( "   ===> in epoch time: ", rdt1, "to", rdt2 )
    print( "       ====> double check: ", epoch2clock(rdt1), "to",  epoch2clock(rdt2))


    if Nforced_stream_length:
        if Nb_min_cnsctv > Nforced_stream_length:
            print('ERROR: `Nb_min_cnsctv` cannot be > `Nforced_stream_length` !'); exit(0)

    
    # Load `distance to coast` data:
    vlon_dist, vlat_dist, xdist = mjt.LoadDist2CoastNC( fdist2coast_nc )

    # Build scan time axis willingly at relative high frequency (dt_bin_sec << dt_buoy_Nmnl)
    NTbin, vTbin, cTbin =   mjt.TimeBins4Scanning( rdt1, rdt2, dt_bin_sec, iverbose=idebug-1 )
    
    # Open, inspect the input file and load raw data:
    Np0, Ns0, vtime0, vykm0, vxkm0, vlat0, vlon0, vBIDs0, vStrm0 = mjt.LoadDataRGPS( cf_in, list_expected_var )
    
    # Masking all point that are before and beyond our period of interest (note: `vtime0` is not masked!):
    Nb, vIDsWP = mjt.KeepDataInterest( rdt1, rdt2, vtime0, vBIDs0, vxkm0, vykm0, vlon0, vlat0,  rmskVal=-99999. )
    # * Nb: number of different buoys that exist for at least 1 record during specified date range aka whole period (WP)
    # * vIDsWP : array(Nb) list (unique) of IDs for these buoys


    mjt.chck4f(cf_npz_strinfo)


    print('\n *** Opening intermediate data into '+cf_npz_strinfo+'!')

    with np.load(cf_npz_strinfo) as data:
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


    # Reminder of what we found with previous script:
    mjt.streamSummaryRGPS( ZNB_ini, ZTc_ini, ZIDs, ZNRc )


    # Loop along streams:
    for js in range(Nstreams):
        
        cs   = str(js)
        vids = np.ma.MaskedArray.compressed( ZIDs[js,:] ) ; # valid IDs for current stream: shrinked, getting rid of masked points
        NvB  = ZNB_ini[js]
        rTc  = ZTc_ini[js]
        cTc  = epoch2clock(rTc)

        if NvB != len(vids): print('ERROR Z1!'); exit(0)
        print('\n *** Having a look at stream #'+cs+' initiated for time bin centered around '+cTc+' !')

        if Nforced_stream_length:
            NCRmax = Nforced_stream_length
        else:        
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
        xxkm = np.zeros((NCRmax,NvB)) - 9999.
        xykm = np.zeros((NCRmax,NvB)) - 9999.
        xlon = np.zeros((NCRmax,NvB)) - 9999.
        xlat = np.zeros((NCRmax,NvB)) - 9999.
        xtim = np.zeros((NCRmax,NvB)) - 9999.      ; # the exact time for each buoy!
        xmsk = np.zeros((NCRmax,NvB) , dtype=int)  ; # the mask for exluding buoys that stick out in time...
        xix0 = np.zeros((NCRmax,NvB) , dtype=int)  ;

        for jb in range(NvB):
            #
            nvr = ZNRc[js,jb] ; # how many successive valid records for this buoy (at least `Nb_min_cnsctv`)
            if Nforced_stream_length:
                nvr = min( nvr, Nforced_stream_length )
            #
            if nvr<Nb_min_cnsctv: print('ERROR Z2!'); exit(0)
            #
            xix0[0:nvr,jb] = ZIX0[js,jb,0:nvr] ; #lolo # all consecutive point position (indices as in `*0` arrays) for thi buoy
            #
            #(idx0_id,) = np.where( vBIDs0 == vids[jb])
            #indv = idx0_id[0:nvr] ; # from record `nvr` onward buoy has been canceled (due to rogue time / expected time)
            indv = xix0[0:nvr,jb].copy() # 
            #
            xxkm[0:nvr,jb] = vxkm0[indv]
            xykm[0:nvr,jb] = vykm0[indv]
            xmsk[0:nvr,jb] = 1
            xlon[0:nvr,jb] = vlon0[indv]
            xlat[0:nvr,jb] = vlat0[indv]
            xtim[0:nvr,jb] = vtime0[indv]


        nBpR[:] = [ np.sum(xmsk[jr,:]) for jr in range(NCRmax) ] ; # How many buoys still present at each record?
        if np.max(nBpR) != NvB: print('ERROR: max(nBpR) != NvB !'); exit(0)
        if idebug>1: print('     +++ num. of boys still present at each record of stream #'+cs+':',nBpR[:])

        # A "mean" time axis along records based on mean accross the buoys remaining at this record...
        ztim = xtim.copy()
        ztim = np.ma.masked_where( xmsk==0, ztim ) ; # otherwize the `mean` in next line would use zeros!!!!
        vtim = np.mean(ztim, axis=1) ; # average on the buoy axis, so `vtim` only dimension is records...
                
        # Nearest interpolation of vtim on the VTbin calendar !!!
        #    => so VT contains the mean date for all buoys at a given record but
        #       corresponding to a value taken from VTbin
        VT = np.zeros( (NCRmax,3), dtype=int )
        i=0
        for rt in vtim:
            idx = np.argmin( np.abs(vTbin[:,0]-rt) )
            if idebug>0: print('      rt =',rt,' => ',epoch2clock(rt),' => nearest of VTbin =',epoch2clock(vTbin[idx,0]))
            VT[i,:] = vTbin[idx,:]
            i=i+1

        # Now, in each record of the stream we should exclude buoys which time position is not inside the expected time bin
        # or is just too far away from the mean of all buoys
        # => if such a buoy is canceld at stream # k, it should also be canceled at following records
        iFU, xmsk, nBpR = mjt.StreamTimeSanityCheck( cs, ztim, VT, xmsk, nBpR, max_t_dev_allowed_in_bin, iverbose=idebug )
        #
        del ztim

        if iFU>0:
            print('old shape =', np.shape(xmsk))                    
            # => we masked some first and/or second buoy records, so we can shrink the arrays accordingly
            (idx_keep0,) , (idx_keep1,) = np.where(xmsk[0,:]==1) , np.where(xmsk[1,:]==1)
            idx_keep = np.unique( np.concatenate([idx_keep0,idx_keep1]) )
            xmsk = xmsk[:,idx_keep]
            xxkm = xxkm[:,idx_keep] 
            xykm = xykm[:,idx_keep] 
            xlon = xlon[:,idx_keep] 
            xlat = xlat[:,idx_keep] 
            xtim = xtim[:,idx_keep] 
            print('new shape =', np.shape(xmsk))


        # Masking:
        xxkm = np.ma.masked_where( xmsk==0, xxkm )
        xykm = np.ma.masked_where( xmsk==0, xykm )
        xlon = np.ma.masked_where( xmsk==0, xlon )
        xlat = np.ma.masked_where( xmsk==0, xlat )
        xtim = np.ma.masked_where( xmsk==0, xtim )
            
        if iplot>0:
            # Stream time evolution on Arctic map:
            kf = mjt.ShowBuoysMap_Trec( vtim, xlon, xlat, pvIDs=[], cnmfig='SELECTION/geo_buoys_RGPS_S'+'%3.3i'%(js),
                                        clock_res='d', NminPnts=Nb_min_buoys )

        # Saving 1 file per stream and per record:
        for jr in range(NCRmax):

            Nbuoys = nBpR[jr] ; # number of buoys alive in current record of this stream
            # Dates based on the center of time bins used:
            idate_binC =         int(VT[jr,0])
            cdate_binC = epoch2clock(VT[jr,0])
            # Dates based on time average accros all buoys at current record of stream
            idate_strM =         int(vtim[jr])
            cdate_strM = epoch2clock(vtim[jr])

            if Nbuoys >= Nb_min_buoys:
                ctr = str.replace( str.replace(cdate_binC[0:16],':','h') ,'-','' ) ; # precision to the minute without ':'
                cf_out = './npz/SELECTION_buoys_RGPS_S'+'%3.3i'%(js)+'_'+ctr+'.npz'

                (indV,) = np.where(xmsk[jr,:]==1) ; # index of valid points
                if len(indV)!= Nbuoys: print('ERROR: rec.',jr,' of stream '+cs+' => len(indV)!= Nbuoys) !', len(indV), Nbuoys); exit(0)

                np.savez_compressed( cf_out, itime=idate_binC, date=cdate_binC, Npoints=Nbuoys, vids=vids[indV],
                                     vtime=xtim[jr,indV], vx=xxkm[jr,indV], vy=xykm[jr,indV], vlon=xlon[jr,indV], vlat=xlat[jr,indV] )

                if iplot>1:
                    # Plot on cartesian coordinates (km):
                    cfpng = './figs/SELECTION/xykm_buoys_RGPS_S'+'%3.3i'%(js)+'_'+ctr+'.png'
                    if jr==0:
                        zrx = [ np.min(xxkm[jr,indV])-100. , np.max(xxkm[jr,indV])+100. ]
                        zry = [ np.min(xykm[jr,indV])-100. , np.max(xykm[jr,indV])+100. ]
                    kg = mjt.ShowTQMesh( xxkm[jr,indV], xykm[jr,indV], ppntIDs=vids[indV], cfig=cfpng, lGeoCoor=False,
                                         zoom=5, rangeX=zrx, rangeY=zry )
            else:
                if idebug>0:
                    print('     ===> NOT saving record #'+str(jr)+' of stream #'+cs+
                          ' (unsufficient n. of buoys alive:',Nbuoys,')')



    ### for js in range(Nstreams)
