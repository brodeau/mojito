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
from scipy import interpolate

from climporn import epoch2clock, clock2epoch
import mojito as mjt


istream = 2 ; # if>0 => restrain the analysis to buoys of stream #istream !

l_time_offset_families = True

cdt_pattern = 'YYYY-MM-DD_00:00:00' ; # pattern for dates

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!

idebug = 0
idebug_interp = 0
iplot = 1          ; # show the result in figures?

list_expected_var = [ 'index', 'lat', 'lon', 'q_flag', 'time' ]

interp_1d = 1 ; # Time interpolation to fixed time axis: 0=> NN ; 1 => linear; 2 => akima

l_drop_coastal   = False ; # get rid of buoys to close to land
MinDistFromLand  = 100  ; # how far from the nearest coast should our buoys be? [km]

l_drop_doublons = True ; # PR: keep the one with the longest record...
rd_tol = 7. # tolerance distance in km to conclude it's the same buoy
NbPass = 2  # number of passes...

l_drop_overlap = False
rhsskm = 7. ; # [km] get rid of buoys (1 of the 2) that are closer to each other than this `rhsskm` km

FillValue = -9999.


dt_nom = 72.2*3600. ; # [s] Nominal `dt` between 2 consecutive records for a buoy => ~ 3 days
pmdt   =  2. *3600. ; # [s] Allowed deviation from `dt_nom` for 2 consecutive records for a buoy => 2 hours



if __name__ == '__main__':

    for cd in ['figs']:
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

    cf_out = 'RGPS_tracking_'+split('_', cdt1)[0]+'_'+split('_', cdt2)[0]+'_lb.nc' ;# netCDF file to generate

    if l_drop_overlap and l_drop_doublons:
        print(' ERROR: you cannot use `l_drop_overlap` and `l_drop_doublons`! Choose one of the two!!!')
        exit(0)

    print('\n *** Date range to restrain data to:')
    print(' ==> '+cdt1+' to '+cdt2 )

    rdt1, rdt2 = clock2epoch(cdt1), clock2epoch(cdt2)
    print( '   ===> in epoch time: ', rdt1, 'to', rdt2 )
    print( '       ====> double check: ', epoch2clock(rdt1), 'to',  epoch2clock(rdt2))

    rdtI = rdt1 + 0.5*dt_bin
    print('\n *** We shall not select buoys that do not make it to at least',epoch2clock(rdtI))

    # Important we want the bins to be centered on the specified dates, so:
    rdtB1 = rdt1 - 0.5*dt_bin

    # Build scan time axis willingly at relative high frequency (dt_bin << dt_buoy_Nmnl)
    Nt, vTbin, cTbin =   mjt.TimeBins4Scanning( rdtB1, rdt2, dt_bin, iverbose=0 )
    if idebug>0:
        for jt in range(Nt):
            print('   --- jt: '+str(jt)+' => ',epoch2clock(vTbin[jt,0]),epoch2clock(vTbin[jt,1]),epoch2clock(vTbin[jt,2]))

    # Now, for proper interpolation right from the start we request any selected buoy to have at least one position
    # under the period [rdtA1:rdt1]:
    rdtA1 = rdt1 - dt_bin
    #rdtA1 = rdt1 - 1.5*dt_bin

    # Open, inspect the input file and load raw data:
    Np0, Ns0, vtime0, vYkm0, vXkm0, vlat0, vlon0, vBIDs0, vStrm0 = mjt.LoadDataRGPS( cf_in, list_expected_var )

    if idebug>0:
        rlat_min, rlat_max = np.min(vlat0), np.max(vlat0)
        print('\n *** Southernmost & Northernmost latitudes in the dataset =>', rlat_min, rlat_max)        
        
    vIDs, idxUNQid = np.unique(vBIDs0, return_index=True )
    vStr           = vStrm0[idxUNQid]
    Nb   = len(vIDs)
    print('\n *** We IDed '+str(Nb)+' (unique) buoys in this file...')
    vmask = np.ones(Nb, dtype='i1')  ; # to mask IDs we cancel

    # First we are going to get rid of all buoys that are not ... our defined period
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    # Stream restriction:
    if istream > 0:
        if istream > Ns0:
            print('ERROR: there are only',Ns0,'streams in the file!!! =>',istream); exit(0)
        #
        (idx_cncl,) = np.where(vStr!=istream)
        vmask[idx_cncl] = 0
        vIDs = np.ma.masked_where( vmask==0, vIDs )
        vIDs = np.ma.MaskedArray.compressed(vIDs)
        vStr = np.ma.masked_where( vmask==0, vStr )
        vStr = np.ma.MaskedArray.compressed(vStr)
        Nb   = len(vIDs)
        vmask = np.ones(Nb, dtype='i1')  ; # to mask IDs we cancel
        print('\n *** UPDATE: based on limiting data to stream '+str(istream)+': '+str(Nb)+' buoys left to follow!')

        
    # For each buoy, getting date of earliest and latest snapshot
    print('\n *** Scanning to retain buoys of interest, based on time range...')

    for jb in range(Nb):
        if jb%1000==0: print('    .... '+str(jb)+' / '+str(Nb)+'...')

        jid = vIDs[jb]

        (idx,) = np.where( vBIDs0==jid )
        vt  = vtime0[idx]
        t1, t2  = np.min(vt), np.max(vt)

        # Get rid of buoys that pop up for the 1st time after rdtB1:
        if t1>rdt1 and vmask[jb]==1:
            vmask[jb] = 0
            if idebug>1: print('   --- excluding buoy #'+str(jb)+' with ID: ',jid,' (pops up after '+epoch2clock(rdt1)+'!)')
        if t2<rdtI and vmask[jb]==1:
            vmask[jb] = 0
            if idebug>1: print('   --- excluding buoy #'+str(jb)+' with ID: ',jid,' (vanishes before '+epoch2clock(rdtI)+'!)')
            
        # Get rid of buoys that do not have at least 1 record point between rdtA1 and rdt1 (important for interpolation onto Vtbin[0,0]!)
        if vmask[jb]==1:
            (idxtp,) = np.where( (vt>rdtA1) & (vt<=rdt1) )
            if len(idxtp) == 0:
                vmask[jb] = 0
                if idebug>1: print('   --- excluding buoy #'+str(jb)+' with ID: ',jid,' (no pos. between:'
                                    +epoch2clock(rdtA1)+' & '+epoch2clock(rdt1)+')')
            
        # Get rid of buoys that do not have a unique record point during the first bin:
        if vmask[jb]==1:
            (idxtp,) = np.where( (vt>vTbin[0,1]) & (vt<=vTbin[0,2]) )
            if len(idxtp) == 0:
                vmask[jb] = 0
                if idebug>1: print('   --- excluding buoy #'+str(jb)+' with ID: ',jid,' (no pos. in 1st time bin:'
                                    +epoch2clock(vTbin[0,1])+' & '+epoch2clock(vTbin[0,2])+')')
            
        if idebug>0 and Nb<=20: print('      * buoy #'+str(jb)+' (id='+str(jid)+': '+epoch2clock(t1)+' ==> '+epoch2clock(t2))


    # Dropping canceled (masked) buoys and updating number of valid buoys:
    vIDs = np.ma.masked_where( vmask==0, vIDs )
    vIDs = np.ma.MaskedArray.compressed(vIDs)
    Nb   = len(vIDs)
    print('\n *** UPDATE: based on date selection, there are '+str(Nb)+' buoys left to follow!')
    

    if l_drop_coastal:
        # For now, only at start time...
        cdata_dir = environ.get('DATA_DIR')
        if cdata_dir==None:
            print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
        fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'
        vlon_dist, vlat_dist, xdist = mjt.LoadDist2CoastNC( fdist2coast_nc ) ; # load `distance to coast` data
        #
        vmask = np.zeros(Nb, dtype='i1') + 1
        print('\n *** Scanning to get rid of buoys too close to land (<'+str(int(MinDistFromLand))+' km)')
        for jb in range(Nb):
            if jb%1000==0: print('    .... '+str(jb)+' / '+str(Nb)+'...')
            jid = vIDs[jb]
            (idx,) = np.where( vBIDs0==jid )
            vt  = vtime0[idx]
            ii = np.argmin( np.abs(vt - rdtB1) ) ; # when we are closest to start time
            zlat, zlon = vlat0[idx[ii]], vlon0[idx[ii]]
            rd_ini = mjt.Dist2Coast( zlon, zlat, vlon_dist, vlat_dist, xdist )
            if rd_ini < MinDistFromLand:
                vmask[jb] = 0
                if idebug>1:
                    print('Buoy '+str(jid)+' too close to land ('+str(int(round(rd_ini,0)))+' km); Lat,Lon ='
                          ,round(zlat,4),',',round(mjt.degE_to_degWE(zlon),4))
            #
        # Dropping canceled (masked) buoys and updating number of valid buoys:
        vIDs = np.ma.masked_where( vmask==0, vIDs      )
        vIDs = np.ma.MaskedArray.compressed(vIDs)
        Nb   = len(vIDs)
        del vmask
        print('\n *** UPDATE: based on minimum distance to nearest coast, there are '+str(Nb)+' buoys left to follow!')
        #
    ### if l_drop_coastal


    
    if l_time_offset_families:
    
        # Analyse at the time calendar of selected buoys
        z1stUsedDate = []
        t_bA =  vTbin[0,1] ; # time at left boundary of 1st target time bin
        ic = 0
        for jid in vIDs:
            ic = ic+1
            (idx,) = np.where( vBIDs0==jid )
            vt   = vtime0[idx] ;            # that's the time axis of this particular buoy [epoch time]
            Np   = len(vt)     ;            # mind that Np varies from 1 buoy to another...
    
            # Inspection of time records for this particular buoy:        
            zdt = vt[1:Np]-vt[0:Np-1]
            if np.any( zdt <= 0.):
                print('ERROR: time not increasing for buoy ID:'+str(jid)+' !!!')
                exit(0)
    
            # Spot the first record of use of the vt series:
            (idxtp,) = np.where( (vt>rdtA1) & (vt<=rdt1) )
            if len(idxtp)==1:
                jt1 = idxtp[0]
            else:
                #print(' AHH: idxtp=',idxtp,' => ',end='')
                #for i in idxtp: print(epoch2clock(vt[i]),'; ',end='')
                #print('')
                idxnrst = np.argmin(abs(vt[idxtp] - t_bA))
                jt1 = idxtp[idxnrst]
            #print('LOLO: jt1 =',jt1,', vt[jt1] =', epoch2clock(vt[jt1]),'\n')
            z1stUsedDate.append(vt[jt1])
    
        z1stUsedDate = np.array(z1stUsedDate,dtype=int)
    
        dt_sens = 600. ; #seconds !
        tfamily = []
        tfamily.append(z1stUsedDate[0])
        for jb in range(Nb):
            lfound = False
            itim = z1stUsedDate[jb]
            #print(epoch2clock(itim))
            for tknown in tfamily:
                if itim>=tknown-dt_sens and itim<tknown+dt_sens:
                    # We know this time... break the loop along tfamily
                    lfound = True
                    break
            if not lfound:
                # It is a time unknown so far:
                print('New time:',epoch2clock(itim))
                tfamily.append(itim)
            
        tfamily = np.array(tfamily,dtype=int)
        ntF = len(tfamily)
        print('We have '+str(ntF)+' different time family identified:')
        for itim in tfamily:
            print(epoch2clock(itim))
    
        ibelongF = np.zeros(Nb, dtype='i1') ; # tells to what family a buoy belongs to...
        for jb in range(Nb):
            itim = z1stUsedDate[jb]
            jf = 0
            for tknown in tfamily:
                if itim>=tknown-dt_sens and itim<tknown+dt_sens:
                    ibelongF[jb] = jf
                    break
                jf = jf+1
    
        vfam = np.zeros(ntF, dtype=int)
        for jf in range(ntF):
            (idxF,) = np.where( ibelongF[:] == jf )
            vfam[jf] = len(idxF)
            print(' * Family #'+str(jf)+' has '+str(vfam[jf])+' buoys!')
    
        kFw = np.argmax( vfam )
        print(' Family '+str(kFw)+' wins! => ', epoch2clock(tfamily[kFw]) )
        if np.sum(vfam) != Nb:
            print('ERRO: `sum(vfam) != Nb` !!!'); exit(0)
    
        # Correction offset for all families:
        vtoffset = np.zeros(ntF, dtype=int)
        vtoffset[:] = tfamily[kFw] - tfamily[:]
        for jf in range(ntF):
            print(' time correction offset to add for family #'+str(jf)+' =>',vtoffset[jf]/3600.,'hours')
            
    ### if l_time_offset_families
    #exit(0)
    
    # Final arrays have 2 dimmensions Nb & Nt (buoy ID & time record)
    xmsk = np.zeros((Nt,Nb), dtype='i1')
    xlon = np.zeros((Nt,Nb)) + FillValue
    xlat = np.zeros((Nt,Nb)) + FillValue
    xXkm = np.zeros((Nt,Nb)) + FillValue
    xYkm = np.zeros((Nt,Nb)) + FillValue

    jb = -1
    for jid in vIDs:
        jb = jb + 1

        (idx,) = np.where( vBIDs0==jid )
        vt   = vtime0[idx] ;            # that's the time axis of this particular buoy [epoch time]
        Np   = len(vt)     ;            # mind that Np varies from 1 buoy to another...

        # Inspection of time records for this particular buoy:        
        zdt = vt[1:Np]-vt[0:Np-1]
        if np.any( zdt <= 0.):
            print('ERROR: time not increasing for buoy ID:'+str(jid)+' !!!')
            exit(0)
        #zmsk = np.zeros(Np-1, dtype='i1') + 1
        #(idxmsk,) = np.where( (zdt>dt_nom+pmdt) | (zdt<dt_nom-pmdt)  )
        #zmsk[idxmsk] = 0
        #print('\nLOLO `dt` in hours:',zdt/3600.) ; # hours
        #print('  => mask =', zmsk)
        #exit(0)

        if l_time_offset_families:
            ifml = ibelongF[jb]
            zadd_offset = vtoffset[ifml]
            #print('LOLO: buoy ID:'+str(jid)+' belongs to family '+str(ifml)+' correction add offset =',zadd_offset/3600.,'hours')
            vt[:] = vt[:] + zadd_offset
        
        ##t_c  = vTbin[0,0] ; # time at center of 1st target time bin
        #t_bA =  vTbin[0,1] ; # time at left boundary of 1st target time bin
        #idxnrst = np.argmin(abs(vt - t_bA)) ; # closest point
        #t_nrst = vt[idxnrst]
        ##print('LOLO: center first time bin and closest point for buoy '+str(jid)+':',epoch2clock(t_c),epoch2clock(t_nrst))
        #print('LOLO: A-bound first time bin and closest point for buoy '+str(jid)+':',epoch2clock(t_bA),epoch2clock(t_nrst))
        ## Check if inside the bin:
        ##if not ( t_nrst>vTbin[0,1] and t_nrst<=vTbin[0,2] ):
        ##    print('ERROR: `t_nrst` not found in first bin!'); exit(0)
        #t_offset = t_bA - t_nrst
        #if idebug>-1: print('   * applying time offset to calendar of buoy '+str(jid)+':', t_offset/3600.,'hours')
        #vt[:] = vt[:] + t_offset


        Nt_b = Nt

        # Does this buoy survive beyond the fixed time range or not??? lilo
        rtend = np.max(vt)
        l_shorten = ( rtend < rdt2 )
        if l_shorten:
            for Nt_b in range(Nt):
                if vTbin[Nt_b,0] > rtend: break
            if idebug>1:
                print('    * buoy :'+str(jid)+' does not make it to '+epoch2clock(rdt2))
                print('      => dies at '+epoch2clock(rtend))
                print('      ==> last index of vTbin to use is: '+str(Nt_b-1)+' => '+epoch2clock(vTbin[Nt_b-1,0]))

        if   interp_1d==0:
            fG = interpolate.interp1d(           vt, vlat0[idx], kind='nearest' )
            fC = interpolate.interp1d(           vt, vYkm0[idx], kind='nearest' )
        elif interp_1d==1:
            fG = interpolate.interp1d(           vt, vlat0[idx], kind='linear' )
            fC = interpolate.interp1d(           vt, vYkm0[idx], kind='linear' )
        elif interp_1d==2:
            fG = interpolate.Akima1DInterpolator(vt, vlat0[idx])
            fC = interpolate.Akima1DInterpolator(vt, vYkm0[idx])
        xlat[0:Nt_b,jb] = fG(vTbin[0:Nt_b,0])
        xYkm[0:Nt_b,jb] = fC(vTbin[0:Nt_b,0])

        if   interp_1d==0:
            fG = interpolate.interp1d(           vt, vlon0[idx], kind='nearest' )
            fC = interpolate.interp1d(           vt, vXkm0[idx], kind='nearest' )
        elif interp_1d==1:
            fG = interpolate.interp1d(           vt, vlon0[idx], kind='linear' )
            fC = interpolate.interp1d(           vt, vXkm0[idx], kind='linear' )
        elif interp_1d==2:
            fG = interpolate.Akima1DInterpolator(vt, vlon0[idx])
            fC = interpolate.Akima1DInterpolator(vt, vXkm0[idx])
        xlon[0:Nt_b,jb] = fG(vTbin[0:Nt_b,0])
        xXkm[0:Nt_b,jb] = fC(vTbin[0:Nt_b,0])

        xmsk[0:Nt_b,jb] = 1

        #if ifml!=kFw: xmsk[0:Nt_b,jb] = 0 ; # LOLO removing all buoys not from winning family!
        
        if idebug_interp>0:
            # Visual control of interpolation:
            mjt.plot_interp_series( jid, 'lat', vt, vTbin[0:Nt_b,:], vlat0[idx], xlat[0:Nt_b,jb] )
            #mjt.plot_interp_series( jid, 'lon', vt, vTbin[0:Nt_b,:], vlon0[idx], xlon[0:Nt_b,jb] )

    if idebug_interp>0: exit(0)
    
    xlat = np.ma.masked_where( xmsk==0, xlat )
    xlon = np.ma.masked_where( xmsk==0, xlon )
    xYkm = np.ma.masked_where( xmsk==0, xYkm )
    xXkm = np.ma.masked_where( xmsk==0, xXkm )


    if l_drop_overlap:
        print('\n *** Applying a spatial downsampling at the scale of '+str(rhsskm)+' km')
        #for jt in range(Nt):
        jt = 0 ; # Only at start !!!!
        #
        (idxn,) = np.where( xmsk[jt,:]==1 )
        zcoor = np.array([ xXkm[jt,:], xYkm[jt,:] ]).T ; # for `gudhi` coord = [X,Y] !
        _, _, idx_keep = mjt.SubSampCloud( rhsskm, zcoor )
        idx_rm = np.setdiff1d( idxn, idx_keep ) ; # keep values of `idxn` that are not in `idx_keep`
        xmsk[jt:,idx_rm] = 0 ; # supressing at this time records and all those following!!!
        #
        Nbn = np.sum(xmsk[0,:])
        print('      => we removed '+str(Nb-Nbn)+' buoys already at first record!')
        (idxK,) = np.where(xmsk[0,:]==1)
        vIDs =   vIDs[idxK]
        xmsk = xmsk[:,idxK]
        xlat = xlat[:,idxK]
        xlon = xlon[:,idxK]
        xYkm = xYkm[:,idxK]
        xXkm = xXkm[:,idxK]
        Nb   = Nbn
        del zcoor, idxn, idx_keep, idx_rm, idxK
        print('\n *** UPDATE: based on "too close to each other" cleaning, there are '+str(Nb)+' buoys left to follow!')



    if l_drop_doublons:
        # Based on first time step only!
        jt = 0
        #
        for jp in range(NbPass):
            print('\n *** Applying initial overlap cleaning at the scale of '+str(rd_tol)+' km')
            zbmask = np.zeros(Nb, dtype='i1') + 1
            zlat,zlon = xlat[jt,:],xlon[jt,:]
            #
            for jb in range(Nb):
                if jb%1000==0: print('       pass #'+str(jp+1)+'.... '+str(jb)+' / '+str(Nb)+'...')
                # All the following work should be done if present buoy has not been cancelled yet
                if zbmask[jb] == 1:
                    jid = vIDs[jb]
                    if idebug>2:
                        print('\n *** "TOO CLOSE" scanning: Buoy #'+str(jb)+'=> ID ='+str(jid))
                        print('    ==> lat,lon = ',zlat[jb],zlon[jb])                
                    vdist = mjt.Haversine( zlat[jb], zlon[jb], zlat, zlon ) ; # build estimates of distances (km) with all other buoy at this same time record:
            
                    # Need to mask itself (distance = 0!)
                    if vdist[jb] != 0.:
                        print(' PROBLEM: distance with yourself should be 0!'); exit(0)
                    vdist[jb] = 9999.
                    if idebug>2:
                        for rd in vdist:
                            if rd<10.: print('   distance = ', rd)
                            
                    rdmin = np.min(vdist)
                    if idebug>2: print('    ==> closest neighbor buoy is at '+str(round(rdmin,2))+' km !')
                    
                    if rdmin < rd_tol:
                        # There is at least 1 buoy too close!
                        # We only deal wit 2 buoys at the time the one we are dealing with and the closest one (otherwize could remove to many)
                        (idx_2c,) = np.where( vdist == rdmin )
                        j2c = idx_2c[0]
                        # Now we have to look at the one of the two that has the longest record:
                        nr1 = np.sum(xmsk[:,jb])
                        nr2 = np.sum(xmsk[:,j2c])
                        if idebug>2: print(' Nb. of valid records for this buoy and the one too close:',nr1,nr2)
                        if nr1<nr2: j2c=jb ; # we should cancel this one not the found one!
                        zbmask[j2c] = 0
                        xmsk[jt:,j2c] = 0

            del zbmask, zlat, zlon
        ### for jp in range(NbPass)

        Nbn = np.sum(xmsk[jt,:])
        print('      => we removed '+str(Nb-Nbn)+' buoys already at first record!')
        (idxK,) = np.where(xmsk[jt,:]==1)
        vIDs =   vIDs[idxK]
        xmsk = xmsk[jt:,idxK]
        xlat = xlat[jt:,idxK]
        xlon = xlon[jt:,idxK]
        xYkm = xYkm[jt:,idxK]
        xXkm = xXkm[jt:,idxK]
        Nb   = Nbn
        del idxK
        
        print('\n *** UPDATE: based on "almost-overlap" cleaning, there are '+str(Nb)+' buoys left to follow! (pass #'+str(jp+1)+')')
    ### if l_drop_doublons

    
    # Masking arrays:
    xlat = np.ma.masked_where( xmsk==0, xlat )
    xlon = np.ma.masked_where( xmsk==0, xlon )
    xYkm = np.ma.masked_where( xmsk==0, xYkm )
    xXkm = np.ma.masked_where( xmsk==0, xXkm )

    # GENERATION OF COMPREHENSIVE NETCDF FILE:
    kk = mjt.ncSaveCloudBuoys( cf_out, vTbin[:,0], vIDs, xYkm, xXkm, xlat, xlon, mask=xmsk, tunits=ctunits_expected, fillVal=FillValue, corigin='RGPS' )

    if iplot>0:

        for jt in range(Nt):
            rt = vTbin[jt,0]
            ct = epoch2clock(rt)
            print(' Ploting for '+ct+'!')
            cfig = './figs/SELECTION/t-interp_buoys_RGPS_'+ct+'.png'
            mjt.ShowBuoysMap( vTbin[jt,0], xlon[jt,:], xlat[jt,:], pvIDs=vIDs, cfig=cfig,
                              ms=5, ralpha=0.5, lShowDate=True, zoom=1., title='RGPS' )
            print('')
