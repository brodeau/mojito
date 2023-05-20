#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

fn_file_dist2coast = 'dist2coast/dist2coast_4deg_North.nc'


known_modes = ['thorough','model','rgps','xlose']


def _check_mode_( caller, mode ):
    if not mode in known_modes:
        print('ERROR [config.'+caller+'()]: unknow mode: '+mode+'!')
        exit(0)
    return 0
    

def initialize( mode='model' ):
    '''
    '''
    from os import environ
    #
    global rc_day2sec, nc_min_buoys_in_batch, nc_forced_batch_length, nc_min_cnsctv, rc_dev_dt_Nmnl, nc_min_buoys_in_batch
    global lc_drop_overlap, rc_Dtol_km
    global rc_Tang_min, rc_Tang_max, rc_Qang_min, rc_Qang_max, rc_dRatio_max
    global lc_accurate_time
    global data_dir, fdist2coast_nc, nc_MinDistFromLand
    
    lk = _check_mode_( 'initialize', mode )

    rc_day2sec = 24*3600

    nc_min_buoys_in_batch = 10 ; # minimum number of buoys for considering a batch a batch!

    nc_MinDistFromLand  = 100. ; # how far from the nearest coast should our buoys be? [km]
    
    nc_forced_batch_length = 2 ; # enforce the length of a batch (each batch will have a maximum of `nc_forced_batch_length` records)
    nc_min_cnsctv = 2          ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)    

    rc_dev_dt_Nmnl = 6*3600    ; # RGPS: max. allowed dev. from the nominal `dt0_RGPS` (~ 3 days) between 2 consecutive records of buoy [s]
    
    lc_drop_overlap = False ; # Because this is done in Quad generation...
    rc_Dtol_km = 6. # tolerance distance in km below which we cancel 1 of the 2 buoys! => for both `l_drop_tooclose` & `lc_drop_overlap`
    #               # Should be same as `rc_d_ss` for 10km as in "updateConfig4Scale()"
    
    lc_accurate_time=True ; # use the exact time at each vertices of the quadrangles when computing deformations

    # Selection of appropriate quadrangles:
    rc_Tang_min =   5. ; # minimum angle tolerable in a triangle [degree]
    rc_Tang_max = 160. ; # maximum angle tolerable in a triangle [degree]
    
    if mode in ['thorough','model']:
        rc_Qang_min =  40.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 140
        rc_dRatio_max = 0.5 ; # value that `max(h1/h2,h2/h1)-1` should not overshoot! h1 being the "height" and "width" of the quadrangle
    elif mode=='rgps':        
        rc_Qang_min =  40.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 140
        rc_dRatio_max = 0.5 ; # value that `max(h1/h2,h2/h1)-1` should not overshoot! h1 being the "height" and "width" of the quadrangle
    elif mode=='xlose':
        rc_Qang_min =  20.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 160.  ; # maximum angle tolerable in a quadrangle [degree]
        rc_dRatio_max = 2. ; # value that `max(h1/h2,h2/h1)-1` should not overshoot! h1 being the "height" and "width" of the quadrangle
        rc_dev_dt_Nmnl = 24*3600 ; # max. allowed dev. from the nominal `dt0_RGPS` (~ 3 days) between 2 consecutive records of buoy [s]

    data_dir = environ.get('DATA_DIR')
    if data_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)

    fdist2coast_nc = data_dir+'/data/'+fn_file_dist2coast

    return 0



def updateConfig4Scale( res_km,  mode='model' ):
    '''
        * res_km: nominal scale for quadrangles we are dealing with [km]
    '''
    import numpy as np
    #
    global ivc_KnownScales, rc_d_ss
    global rc_Qarea_min, rc_Qarea_max, rc_tolQuadA, rc_t_dev_cancel
    global rc_div_min, rc_shr_min, rc_tot_min, rc_div_max, rc_shr_max, rc_tot_max
    global rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig
    
    irk = int(res_km)

    lk = _check_mode_( 'updateConfig4Scale', mode )
    
    rc_tolQuadA = 3. * res_km/20. ; # +- tolerance in [km] on the MEAN scale of quadrangles in a batch
    #                               #    to accept a given scale. Ex: average scale of quadrangle = 15.9 km is accepted for 15 km !!

    rc_d_ss = 6. ; # default radius for subsampling the cloud of points [km]

    # Extremas for figures only:
    rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig =  1., 1., 1.
    # Extremas for deformations, in [day^-1]:
    rc_div_max, rc_shr_max, rc_tot_max = 1., 1., 1.


    if mode in ['rgps']:
        # Because noisy at tiny deformation!!!        
        rc_div_min, rc_shr_min, rc_tot_min = 8.e-5, 8.e-5, 8.e-5    ; # based on scatter plot of scaling        
    else:
        rc_div_min, rc_shr_min, rc_tot_min = 1.e-15, 1.e-15, 1.e-15 ; # for deformation generation
    
    if   irk==10:
        rc_d_ss = 6.
        rc_tolQuadA = 2
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig =  0.1, 0.05, 0.05
        if mode in ['thorough','model']:            
            #rc_Qarea_min, rc_Qarea_max = 9.875*9.875, 10.125*10.125
            rc_Qarea_min, rc_Qarea_max = 9.75*9.75, 10.25*10.25
        elif mode=='rgps':
            # When lose, average scale is 10.33 km with a standard deviation of ~1km!
            rc_Qarea_min, rc_Qarea_max = 9.83*9.83, 10.83*10.83 ; # 10.33 +- 0.5
        elif mode=='xlose':
            rc_tolQuadA = 10
            rc_Qarea_min, rc_Qarea_max = 5*5, 15*15
        else:
            rc_Qarea_min, rc_Qarea_max = 8.*8., 12.*12.
        #
    elif irk==20:
        rc_d_ss = 15.
        rc_tolQuadA = 4
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig = 5.e-1, 5.e-1, 5.e-1
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max = 19.5*19.5, 20.5*20.5
        elif mode=='rgps':
            rc_Qarea_min, rc_Qarea_max = 19.*19., 21.*21.
            #rc_div_max, rc_shr_max, rc_tot_max = 0.7, 0.7, 0.7 ; # day^-1 ; => supresses irrealistically large values
        else:
            rc_Qarea_min, rc_Qarea_max = 18.*18., 22.*22.
        #
    elif irk==40:
        rc_d_ss = 34.5
        rc_tolQuadA = 8
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig = 2.e-1, 2.e-1, 2.e-1
        if mode=='model':            
            #rc_Qarea_min, rc_Qarea_max = 39.5*39.5, 40.5*40.5
            rc_Qarea_min, rc_Qarea_max = 39.*39., 41.*41.
        elif mode=='rgps':
            rc_Qarea_min, rc_Qarea_max = 38.*38., 42.*42.
            #rc_div_max, rc_shr_max, rc_tot_max = 0.4, 0.4, 0.4 ; # day^-1 ; => supresses irrealistically large values
        else:
            rc_Qarea_min, rc_Qarea_max = 35*35, 45*45
        #
    elif irk==80:
        rc_d_ss = 74.75
        rc_tolQuadA = 16
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig = 1.e-1, 1.e-1, 1.e-1
        if mode=='model':            
            rc_Qarea_min, rc_Qarea_max = 78.*78., 82.*82.
        elif mode=='rgps':
            rc_Qarea_min, rc_Qarea_max = 76.*76., 84.*84.
            #rc_div_max, rc_shr_max, rc_tot_max = 0.2, 0.2, 0.2 ; # day^-1 ; => supresses irrealistically large values
        else:
            rc_Qarea_min, rc_Qarea_max = 70*70, 90*90
        #
    elif irk==160:
        rc_d_ss = 156.
        rc_tolQuadA = 32
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig = 7.5e-2, 7.5e-2, 7.5e-2
        if mode=='model':            
            rc_Qarea_min, rc_Qarea_max = 156*156, 164*164
        elif mode=='rgps':
            rc_Qarea_min, rc_Qarea_max = 152.*152., 168.*168.
            #rc_div_max, rc_shr_max, rc_tot_max = 0.1, 0.1, 0.1 ; # day^-1 ; => supresses irrealistically large values
        else:
            rc_Qarea_min, rc_Qarea_max = 140*140, 180*180
        #
    elif irk==320:
        rc_d_ss = 315.6
        rc_tolQuadA = 64
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig = 5.e-2, 5.e-2, 5.e-2
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max =  312*312, 328*328
        elif mode=='rgps':            
            rc_Qarea_min, rc_Qarea_max = 304.*304., 336.*336
            #rc_div_max, rc_shr_max, rc_tot_max = 6.e-2, 6.e-2, 6.e-2 ; # day^-1 ; => supresses irrealistically large values
        else:
            rc_Qarea_min, rc_Qarea_max =  280*280, 360*360            
        #
    elif irk==640:
        rc_d_ss = 636.
        rc_tolQuadA = 128
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig = 2.5e-2, 2.5e-2, 2.5e-2
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max =  624*624, 656*656
        elif mode=='rgps':            
            rc_Qarea_min, rc_Qarea_max = 608.*608., 672.*672.
            #rc_div_max, rc_shr_max, rc_tot_max = 4.e-2, 4.e-2, 4.e-2 ; # day^-1 ; => supresses irrealistically large values
        else:
            rc_Qarea_min, rc_Qarea_max =  480*480, 800*800
        #
    else:
        print('ERROR [updateConfig4Scale()]: scale "'+str(irk)+' km" is unknown!'); sys.exit(0)

    
    print(' *** [updateConfig4Scale](): upper and lower bound for scale "'+str(res_km)+' km":', np.sqrt([rc_Qarea_min, rc_Qarea_max]))
    print('                             => rc_Qarea_min, Qarea_Nom, rc_Qarea_max =',rc_Qarea_min, res_km*res_km, rc_Qarea_max,'km^2')

    #############################

    # Only for RGPS data:
    # *******************
    # When not at the nominal scale, we can adapt `rc_t_dev_cancel` to the scale we are dealing with:
    #    a quadrangle can involve points that are too distant from one another in terms of time
    #    => we disregard any quadrangles which standard deviation of the time of the 4 positions
    #       excess `t_dev_cancel` seconds !
    rc_t_dev_cancel = 60

    if mode == 'rgps':
        if   res_km>=70. and res_km<150.:
            rc_t_dev_cancel =  3.*3600
        elif   res_km>=150. and res_km<300.:
            rc_t_dev_cancel =  3600   ; #lolo? For now, doesn't make a big diff appart from less def values...
            #rc_t_dev_cancel =  6*3600   ; #ok
            #rc_t_dev_cancel =  12*3600   ; #lolo?
            #rc_t_dev_cancel =  24*3600   ; #lolo?
        elif res_km>=300. and res_km<600.:
            #rc_t_dev_cancel =  12*3600
            #rc_t_dev_cancel =  18*3600
            rc_t_dev_cancel =  24*3600 ; #ok! seems like the acceptable shit! 20230520
        elif res_km>=600.:
            #rc_t_dev_cancel = 6*3600
            #rc_t_dev_cancel = 9*3600
            #rc_t_dev_cancel = 12*3600
            rc_t_dev_cancel = 24*3600 ; # ok
        #
    if rc_t_dev_cancel > 60:
        print('\n *** `rc_t_dev_cancel` updated to ',rc_t_dev_cancel/3600,'hours!')

    return 0
