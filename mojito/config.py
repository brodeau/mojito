#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

fn_file_dist2coast = 'dist2coast/dist2coast_4deg_North.nc'

def initialize():

    global rc_day2sec, nc_min_buoys_in_batch
    global rc_Tang_min, rc_Tang_max, rc_Qang_min, rc_Qang_max, rc_dRatio_max
    global lc_accurate_time, rc_div_max, rc_shr_max, rc_tot_max
    global data_dir, fdist2coast_nc, nc_MinDistFromLand

    rc_day2sec = 24*3600

    nc_min_buoys_in_batch = 10 ; # minimum number of buoys for considering a batch a batch!

    nc_MinDistFromLand  = 100. ; # how far from the nearest coast should our buoys be? [km]


    # Selection of appropriate quadrangles:
    rc_Tang_min =   5. ; # minimum angle tolerable in a triangle [degree]
    rc_Tang_max = 160. ; # maximum angle tolerable in a triangle [degree]
    #
    rc_Qang_min =  40.  ; # minimum angle tolerable in a quadrangle [degree]
    rc_Qang_max = 150.  ; # maximum angle tolerable in a quadrangle [degree]
    rc_dRatio_max = 0.5 ; # value that `max(h1/h2,h2/h1)-1` should not overshoot! h1 being the "height" and "width" of the quadrangle


    lc_accurate_time=True ; # use the exact time at each vertices of the quadrangles when computing deformations

    # For figures, in days^-1:
    rc_div_max = 0.1
    rc_shr_max = 0.1
    rc_tot_max = 0.1        

    data_dir = environ.get('DATA_DIR')
    if cdata_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
        fdist2coast_nc = cdata_dir+'/data/'+fn_file_dist2coast

    return 0




def updateConfig4Scale( res_km, binDt ):
    '''
        * res_km: nominal scale for quadrangles we are dealing with [km]
        * binDt : width of time bins we have used for selection period (6h, 1 day, 3 days?) [s]
    '''
    import numpy as np
    #
    global ivc_KnownScales
    global rc_Qarea_min, rc_Qarea_max, rc_tolQuadA, rc_t_dev_cancel
    

    rdev_scale = 0.1 ; # how much can we deviate from the specified scale to accept or reject a quadrangle
    #                  # =>  (reskm*(1-rdev_scale))**2  <  Quadrangles_area < (reskm*(1+rdev_scale))**2
    
    #zA  = res_km**2
    #dR2 = abs( zA - (res_km*(1+rdev_scale))**2 )
    #rc_Qarea_min, rc_Qarea_max = zA-dR2, zA+dR2
    # Scales to retain Quads:
    vscales = np.array( [ 5*2**i  for i in range(10) ], dtype=int )
    vUb = (vscales[0:-1]+vscales[1:])/2.
    vLb = np.zeros(len(vUb))
    vLb[1:] = vUb[0:-1]
    #
    ivc_KnownScales = vscales[1:-1]
    del vscales
    vLb = vLb[1:]
    vUb = vUb[1:]
    Ns = len(ivc_KnownScales)
    if len(vLb)!=Ns or len(vUb)!=Ns:
        print('ERROR [updateConfig4Scale]: `len(vLb)!=Ns or len(vUb)!=Ns`')
    vtolerc = np.array( [ 0.1*2**i  for i in range(Ns) ] )
    #print('LOLO: toler. =', vtolerc)
    vLb, vUb = vLb + vtolerc , vUb - vtolerc
    #for ks in range(Ns-1):
    #    print('* Scale =',ivc_KnownScales[ks],'  => lower and upper bounds =', vLb[ks], vUb[ks])
    if not int(res_km) in ivc_KnownScales:
        print('ERROR [updateConfig4Scale]: scale "'+str(res_km)+' km" is unknown!'); sys.exit(0)
    ([js],) = np.where(ivc_KnownScales==int(res_km))
    rc_Qarea_min, rc_Qarea_max = vLb[js]*vLb[js], vUb[js]*vUb[js]    
    print(' *** [updateConfig4Scale](): upper and lower bound for scale "'+str(res_km)+' km":', vLb[js],vUb[js])
    print('                             => rc_Qarea_min, Qarea_Nom, rc_Qarea_max =',rc_Qarea_min, res_km*res_km, rc_Qarea_max,'km^2')
        


    #############################
    rc_tolQuadA = 3. * res_km/20. ; # +- tolerance in [km] on the MEAN scale of quadrangles in a batch
    #                               #    to accept a given scale. Ex: average scale of quadrangle = 15.9 km is accepted for 15 km !!
    if res_km>35. and res_km<45.:
        rc_tolQuadA = 5.
    if res_km>70. and res_km<300:
        rc_tolQuadA = 15.
    if res_km>=300.:
        rc_tolQuadA = 50.
    if res_km>=600.:
        rc_tolQuadA = 100.

    # When not at the nominal scale, we can adapt `rc_t_dev_cancel` to the scale we are dealing with:
    #    a quadrangle can involve points that are too distant from one another in terms of time
    #    => we disregard any quadrangles which standard deviation of the time of the 4 positions
    #       excess `t_dev_cancel` seconds !
    rc_t_dev_cancel = 60
    if binDt>6*3600:
        # `rc_t_dev_cancel` remains at 60s when the selection bin with is 6 hours or below
        if res_km>=100.:
            rc_t_dev_cancel =  3600
        if res_km>=150.:
            rc_t_dev_cancel = 1.5*3600
        if res_km>=300.:
            rc_t_dev_cancel = 3*3600
        if res_km>=600.:
            rc_t_dev_cancel = 6*3600
            print('\n *** `rc_t_dev_cancel` updated to ',rc_t_dev_cancel/3600,'hours!')

    exit(0)

            
