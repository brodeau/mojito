#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

fn_file_dist2coast = 'dist2coast/dist2coast_4deg_North.nc'


def initialize( mode='thorough' ):
    '''
    '''
    from os import environ
    #
    global rc_day2sec, nc_min_buoys_in_batch, nc_forced_batch_length, nc_min_cnsctv, rc_dev_dt_Nmnl, nc_min_buoys_in_batch
    global lc_drop_overlap, rc_Dtol_km
    global rc_Tang_min, rc_Tang_max, rc_Qang_min, rc_Qang_max, rc_dRatio_max
    global lc_accurate_time
    global rc_div_min, rc_shr_min, rc_tot_min, rc_div_max, rc_shr_max, rc_tot_max
    global data_dir, fdist2coast_nc, nc_MinDistFromLand

    rc_day2sec = 24*3600

    nc_min_buoys_in_batch = 10 ; # minimum number of buoys for considering a batch a batch!

    nc_MinDistFromLand  = 100. ; # how far from the nearest coast should our buoys be? [km]
    
    nc_forced_batch_length = 2 ; # enforce the length of a batch (each batch will have a maximum of `nc_forced_batch_length` records)
    nc_min_cnsctv = 2        ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)    

    rc_dev_dt_Nmnl = 6*3600  ; # maximum allowed deviation from the nominal `dt0_RGPS` (~ 3 days) between 2 consecutive records of buoy [s]
    
    lc_drop_overlap = True
    rc_Dtol_km = 5.5 # tolerance distance in km below which we decide to cancel one of the 2 buoys! => for both `l_drop_tooclose` & `lc_drop_overlap`

    # Selection of appropriate quadrangles:
    rc_Tang_min =   5. ; # minimum angle tolerable in a triangle [degree]
    rc_Tang_max = 160. ; # maximum angle tolerable in a triangle [degree]
    #
    if mode in ['thorough','model']:
        rc_Qang_min =  40.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 150.  ; # maximum angle tolerable in a quadrangle [degree]
        rc_dRatio_max = 0.5 ; # value that `max(h1/h2,h2/h1)-1` should not overshoot! h1 being the "height" and "width" of the quadrangle
    elif mode=='xlose':
        rc_Qang_min =  30.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 160.  ; # maximum angle tolerable in a quadrangle [degree]
        rc_dRatio_max = 3. ; # value that `max(h1/h2,h2/h1)-1` should not overshoot! h1 being the "height" and "width" of the quadrangle
    else:
        print('ERROR [initialize()]: unknow mode: '+mode+'!'); exit(0)

    lc_accurate_time=True ; # use the exact time at each vertices of the quadrangles when computing deformations
    #lc_accurate_time=False ; # use the exact time at each vertices of the quadrangles when computing deformations

    rc_div_min, rc_shr_min, rc_tot_min = 0.8e-4, 0.8e-4, 0.8e-4   ; # for scaling (not PDFs)
    rc_div_max, rc_shr_max, rc_tot_max = 0.1, 0.1, 0.1            ; # for figures

    data_dir = environ.get('DATA_DIR')
    if data_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)

    fdist2coast_nc = data_dir+'/data/'+fn_file_dist2coast

    return 0




def updateConfig4Scale( res_km,  mode='thorough' ):
    '''
        * res_km: nominal scale for quadrangles we are dealing with [km]
    '''
    import numpy as np
    #
    global ivc_KnownScales, rc_d_ss
    global rc_Qarea_min, rc_Qarea_max, rc_tolQuadA, rc_t_dev_cancel

    irk = int(res_km)

    rc_tolQuadA = 3. * res_km/20. ; # +- tolerance in [km] on the MEAN scale of quadrangles in a batch
    #                               #    to accept a given scale. Ex: average scale of quadrangle = 15.9 km is accepted for 15 km !!

    rc_d_ss = 8 ; # radius for subsampling the cloud of points [km]
    
    
    min_div, min_shr, min_tot = 0.003, 0.003, 0.003 ; # day^-1 ; RGPS is noisy around 0! We do not want have the zero on the PDF...
    
    if   irk==10:
        if mode in ['thorough','model']:
            rc_tolQuadA = 1
            rc_Qarea_min, rc_Qarea_max = 9*9, 11*11
        elif mode=='xlose':
            rc_tolQuadA = 5
            rc_Qarea_min, rc_Qarea_max = 2*2, 20*20
        else:
            print('ERROR [updateConfig4Scale()]: unknow mode: '+mode+'!'); exit(0)

    elif irk==20:
        rc_d_ss = 14.6
        rc_tolQuadA = 3
        rc_Qarea_min, rc_Qarea_max = 18*18, 22*22
        min_div, min_shr, min_tot = 0.001, 0.001, 0.001
        #
    elif irk==40:
        rc_d_ss = 34.5
        rc_tolQuadA = 7
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max = 38*38, 42*42
        else:
            rc_Qarea_min, rc_Qarea_max = 35*35, 45*45
        #
    elif irk==80:
        rc_d_ss = 74.75
        rc_tolQuadA = 15
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max = 78*78, 82*82
        else:
            rc_Qarea_min, rc_Qarea_max = 70*70, 90*90
        #
    elif irk==160:
        rc_d_ss = 156.
        rc_tolQuadA = 30
        if mode=='model':        
            rc_Qarea_min, rc_Qarea_max = 155*155, 165*165
        else:
            rc_Qarea_min, rc_Qarea_max = 140*140, 180*180
        #
    elif irk==320:
        rc_d_ss = 315.6
        rc_tolQuadA = 100
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max =  315*315, 325*325
        else:
            ##rc_Qarea_min, rc_Qarea_max =  260*260, 380*380
            rc_Qarea_min, rc_Qarea_max =  280*280, 360*360
            #rc_Qarea_min, rc_Qarea_max =  310*310, 330*330
        #
    elif irk==640:
        rc_d_ss = 636.
        rc_tolQuadA = 300
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max =  635*635, 645*645
        else:
            ##rc_Qarea_min, rc_Qarea_max =  600*600, 680*680
            ##rc_Qarea_min, rc_Qarea_max =  520*520, 760*760
            ##rc_Qarea_min, rc_Qarea_max =  500*500, 780*780
            rc_Qarea_min, rc_Qarea_max =  480*480, 800*800
        #
    else:
        print('ERROR [updateConfig4Scale]: scale "'+str(irk)+' km" is unknown!'); sys.exit(0)

        
  

        
    #vscales = np.array( [ 5*2**i  for i in range(10) ], dtype=int )
    #vUb = (vscales[0:-1]+vscales[1:])/2.
    #vLb = np.zeros(len(vUb))
    #vLb[1:] = vUb[0:-1]
    #ivc_KnownScales = vscales[1:-1]
    #del vscales
    #vLb = vLb[1:]
    #vUb = vUb[1:]
    #Ns = len(ivc_KnownScales)
    #if len(vLb)!=Ns or len(vUb)!=Ns:
    #    print('ERROR [updateConfig4Scale]: `len(vLb)!=Ns or len(vUb)!=Ns`')
    #vtolerc = np.array( [ 0.1*2**i  for i in range(Ns) ] ) ; #print('LOLO: toler. =', vtolerc)
    #vLb, vUb = vLb + vtolerc , vUb - vtolerc
    ##for ks in range(Ns-1):
    ##    print('* Scale =',ivc_KnownScales[ks],'  => lower and upper bounds =', vLb[ks], vUb[ks])
    #if not int(res_km) in ivc_KnownScales:
    #    print('ERROR [updateConfig4Scale]: scale "'+str(res_km)+' km" is unknown!'); sys.exit(0)
    #([js],) = np.where(ivc_KnownScales==int(res_km))
    #rc_Qarea_min, rc_Qarea_max = vLb[js]*vLb[js], vUb[js]*vUb[js]

    print(' *** [updateConfig4Scale](): upper and lower bound for scale "'+str(res_km)+' km":', np.sqrt([rc_Qarea_min, rc_Qarea_max]))
    print('                             => rc_Qarea_min, Qarea_Nom, rc_Qarea_max =',rc_Qarea_min, res_km*res_km, rc_Qarea_max,'km^2')

    #############################

    # When not at the nominal scale, we can adapt `rc_t_dev_cancel` to the scale we are dealing with:
    #    a quadrangle can involve points that are too distant from one another in terms of time
    #    => we disregard any quadrangles which standard deviation of the time of the 4 positions
    #       excess `t_dev_cancel` seconds !
    rc_t_dev_cancel = 60
    if res_km>=300.:
        rc_t_dev_cancel = 12*3600
        if res_km>=600.:
            rc_t_dev_cancel = 24*3600
            #rc_t_dev_cancel = 18*3600
        print('\n *** `rc_t_dev_cancel` updated to ',rc_t_dev_cancel/3600,'hours!')

    return 0
