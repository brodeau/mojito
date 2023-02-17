l_drop_doublons = False ; # PR: keep the one with the longest record...
rd_tol = 5. # tolerance distance in km to conclude it's the same buoy

    # Removing doublons:
    if l_drop_doublons:
        print(' #FIXME: l_drop_doublons !!!'); exit(0)
        
        # Use x,y !!! not lon,lat likebelow !!!
        from .tracking import Haversine

        vid_mask = np.zeros(Nb, dtype='i1')
        vid_mask[:] = 1
        
        # Based on first time step:
        jt   = 0
        vlat = xlat[jt,:]
        vlon = xlon[jt,:]
        #
        for jb in range(Nb):        
            # All the following work should be done if present buoy has not been cancelled yet
            if vid_mask[jb] == 1:
                
                jid = vIDs[jb]
                if idebug>0: print('\n *** "TOO CLOSE" scanning: Buoy #'+str(jb)+'=> ID ='+str(jid))
        
                rlat = vlat[jb]
                rlon = vlon[jb]
        
                if idebug>0: print('    ==> lat,lon = ',rlat,rlon)
                
                # Building array of distances (km) with all other buoy at this same time record:
                vdist = Haversine( rlat, rlon, vlat, vlon )
        
                # Need to mask itself (distance = 0!)
                if vdist[jb] != 0.:
                    print(' PROBLEM: distance with yourself should be 0!')
                vdist[jb] = 9999.
        
                if idebug>1:
                    for rd in vdist:
                        if rd<10.: print('   distance = ', rd)
        
                rdmin = np.min(vdist)
                if idebug>0: print('    ==> closest other buoy is at '+str(round(rdmin,2))+' km !')
                
                if rdmin < rd_tol:
                    # There is at least 1 buoy too close!
                    #  => populating the buoys closer than `rd_tol` km:
                    (vidx_close,) = np.where( vdist < rd_tol )
        
                    if idebug>0:
                        print('    ==> list of buoys closer than '+str(rd_tol)+' km:')
                        for jb in vidx_close:
                            print(vIDs[ii])
    
                    vid_mask[(vidx_close,)] = 0
        
        # Alright, all IDs spotted via vid_mask should be cancelled in the 2D fields, at all time steps
        (idmsk,) = np.where(vid_mask==0)
        xmsk[:,idmsk] = 0
        vIDs = np.ma.masked_where( vid_mask==0, vIDs )
        xlat = np.ma.masked_where( xmsk==0, xlat )
        xlon = np.ma.masked_where( xmsk==0, xlon )
        
        Nok = np.sum(xmsk[0,:])
        print('\n *** UPDATE: after "too close" canceling, there are '+str(Nok)+' buoys left at t0!')

        # Need to compress arrays, didn't find an elegant way to do it so here is the following abomination:
        xmsk2 = np.zeros((Nt,Nok), dtype=int)
        xlon2 = np.zeros((Nt,Nok))
        xlat2 = np.zeros((Nt,Nok))
        ic = 0
        for jb in range(Nb):
            if xmsk[0,jb] == 1:
                xmsk2[:,ic] = xmsk[:,jb]
                xlon2[:,ic] = xlon[:,jb]
                xlat2[:,ic] = xlat[:,jb]
                ic = ic+1
        vIDs = np.ma.MaskedArray.compressed(vIDs)
        Nb = Nok
        del xmsk, xlat, xlon
        xmsk = xmsk2
        xlat = xlat2
        xlon = xlon2
        del xmsk2, xlat2, xlon2        
    ### if l_drop_doublons

    
