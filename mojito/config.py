#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

fn_file_dist2coast = 'dist2coast/dist2coast_4deg_North.nc'


known_modes = ['thorough','model','rgps','xlose','rgps_map','rgps_track']


def controlModeName( caller, mode ):
    if not mode in known_modes:
        print('ERROR [config.controlModeName::'+caller+']: unknow mode: '+mode+'!')
        exit(0)
    return 0
    

def initialize( mode='model' ):
    '''
    '''
    from os import environ
    #
    global rc_day2sec, nc_min_buoys_in_batch, nc_forced_batch_length, nc_min_cnsctv, rc_dev_dt_Nmnl, nc_min_buoys_in_batch
    global lc_drop_overlap, rc_Dtol_km
    global rc_Tang_min, rc_Tang_max
    global lc_accurate_time
    global data_dir, fdist2coast_nc, nc_MinDistFromLand
    
    lk = controlModeName( 'initialize', mode )

    rc_day2sec = 24*3600

    nc_min_buoys_in_batch = 10 ; # minimum number of buoys for considering a batch a batch!

    nc_MinDistFromLand  = 100. ; # how far from the nearest coast should our buoys be? [km]
    
    nc_forced_batch_length = 2 ; # enforce the length of a batch (each batch will have a maximum of `nc_forced_batch_length` records)
    nc_min_cnsctv = 2          ; # minimum number of consecutive buoy positions to store (>=2, because we need to do a d/dt)    

    rc_dev_dt_Nmnl = 6*3600    ; # RGPS: max. allowed dev. from the nominal `dt0_RGPS` (~ 3 days) between 2 consecutive records of buoy [s]
    if mode in ['xlose']:
        rc_dev_dt_Nmnl = 24*3600 ; # max. allowed dev. from the nominal `dt0_RGPS` (~ 3 days) between 2 consecutive records of buoy [s]
    
    lc_drop_overlap = False ; # Because this is done in Quad generation...
    rc_Dtol_km = 5. # (if lc_drop_overlap) tolerance distance in km below which we cancel 1 of the 2 buoys! => for both `l_drop_tooclose` & `lc_drop_overlap`
    #               # Should be same as `rc_d_ss` for 10km as in "updateConfig4Scale()"

    lc_accurate_time=True ; # use the exact time at each vertices of the quadrangles when computing deformations

    # Selection of appropriate quadrangles:
    rc_Tang_min =   5. ; # minimum angle tolerable in a triangle [degree]
    rc_Tang_max = 160. ; # maximum angle tolerable in a triangle [degree]
    

    data_dir = environ.get('DATA_DIR')
    if data_dir==None:
        print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)

    fdist2coast_nc = data_dir+'/data/'+fn_file_dist2coast

    return 0



def updateConfig4Scale( res_km,  mode='model', ltalk=True ):
    '''
        * res_km: nominal scale for quadrangles we are dealing with [km]
    '''
    import numpy as np
    #
    global ivc_KnownScales, rc_d_ss, rcFracOverlapOK
    global rc_Qang_min, rc_Qang_max, rc_dRatio_max
    global rc_Qarea_min, rc_Qarea_max, rc_maxDevMeanAreaQuads, rc_t_dev_cancel
    global rc_div_min, rc_shr_min, rc_tot_min, rc_div_max, rc_shr_max, rc_tot_max
    global rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig
    global rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf, rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf
    
    irk = int(res_km)

    lk = controlModeName( 'updateConfig4Scale', mode )
    
    rc_maxDevMeanAreaQuads = 3. * res_km/20. ; # +- tolerance in [km] on the MEAN scale of quadrangles in a batch
    #                               #    to accept a given scale. Ex: average scale of quadrangle = 15.9 km is accepted for 15 km !!

    rc_d_ss = 6. ; # default radius for subsampling the cloud of points [km]


    rcFracOverlapOK = 0.05 ; # when multi-realisation of coarsening (`rd_ss`) tells how much the two quadrangles from 2 different realisations can overlap (area),
    #                        # as a fraction of the nominal areas (=`reskm*reskm`)
    
    # Extremas for deformations, in [day^-1]:
    rc_div_max, rc_shr_max, rc_tot_max = 1., 1., 1.
    if mode in ['rgps']:
        rc_div_min, rc_shr_min, rc_tot_min = 8.e-5, 8.e-5, 8.e-5    ; # # Because noisy at tiny deformation!!! Based on scatter plot of scaling        
    else:
        rc_div_min, rc_shr_min, rc_tot_min = 1.e-15, 1.e-15, 1.e-15 ; # for deformation generation

    # Extremas for deformations for PDFs only:
    rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.5, 0.5, 0.5
    rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf =  8.e-5, 8.e-5, 8.e-5
    
    # Extremas for deformations for figures only:
    rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig =  1., 1., 1.,0.1




    if mode in ['thorough','model']:
        rc_Qang_min =  40.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 140.
        rc_dRatio_max = 0.5 ; # value that `max(h1/h2,h2/h1)-1` should not exceed! h1 being the "height" and "width" of the quadrangle
        #
    elif mode in ['rgps','rgps_track']:        
        rc_Qang_min =  40.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 140.
        rc_dRatio_max = 0.5 ; # value that `max(h1/h2,h2/h1)-1` should not exceed! h1 being the "height" and "width" of the quadrangle
        #
        if irk in [640]:
            rc_Qang_min =  20. ; # minimum angle tolerable in a quadrangle [degree]
        #
    elif mode in ['rgps_map']:        
        rc_Qang_min =  20.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 160.
        rc_dRatio_max = 0.9 ; # value that `max(h1/h2,h2/h1)-1` should not exceed! h1 being the "height" and "width" of the quadrangle
    elif mode in ['xlose']:
        rc_Qang_min =  20.  ; # minimum angle tolerable in a quadrangle [degree]
        rc_Qang_max = 160.  ; # maximum angle tolerable in a quadrangle [degree]
        rc_dRatio_max = 1. ; # value that `max(h1/h2,h2/h1)-1` should not exceed! h1 being the "height" and "width" of the quadrangle








    
    if mode=='rgps_map' and irk>10:
        print('ERROR [updateConfig4Scale]: "mode=rgps_map" is only meant to be used at 10km !!!')
        exit(0)
    
    if   irk==10:
        rc_d_ss = 5.5
        rcFracOverlapOK = 1.e-9
        rc_maxDevMeanAreaQuads = 2        
        rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf = 0.003, 0.003, 0.003
        rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.2, 0.2, 0.2
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig =  0.12, 0.12, 0.12, 0.05
        if mode in ['thorough','model']:            
            rc_Qarea_min, rc_Qarea_max = 9.75*9.75, 10.25*10.25
        elif mode=='rgps_map':
            rc_maxDevMeanAreaQuads = 5
            rc_Qarea_min, rc_Qarea_max = 6.*6.,14.*14.
        elif mode=='xlose':
            rc_maxDevMeanAreaQuads = 10
            rc_Qarea_min, rc_Qarea_max = 5*5, 15*15
        else:
            rc_Qarea_min, rc_Qarea_max = 8.*8., 12.*12.
        #
    elif irk==20:
        rc_d_ss = 15.
        rcFracOverlapOK = 1.e-9
        rc_maxDevMeanAreaQuads = 4
        rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf = 0.001, 0.001, 0.001
        rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.5, 0.5, 0.5
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig = 1.e-1, 1.e-1, 1.e-1, 0.02
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max = 19.5*19.5, 20.5*20.5
        else:
            rc_Qarea_min, rc_Qarea_max = 18.*18., 22.*22.
        #
    elif irk==40:
        rc_d_ss = 34.5
        rc_maxDevMeanAreaQuads = 8
        rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf = 5.e-4, 5.e-4, 5.e-4
        rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.2, 0.2, 0.2
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig = 0.5e-1, 0.5e-1, 0.5e-1, 0.01
        if mode=='model':            
            #rc_Qarea_min, rc_Qarea_max = 39.5*39.5, 40.5*40.5
            rc_Qarea_min, rc_Qarea_max = 39.*39., 41.*41.
        else:
            rc_Qarea_min, rc_Qarea_max = 35*35, 45*45
        #
    elif irk==80:
        rc_d_ss = 74.75
        rcFracOverlapOK = 0.05
        rc_maxDevMeanAreaQuads = 16
        rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf = 1.e-4, 1.e-4, 1.e-4
        rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.15, 0.15, 0.15
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig = 1.e-1, 1.e-1, 1.e-1, 0.025
        if mode=='model':            
            rc_Qarea_min, rc_Qarea_max = 78.*78., 82.*82.
        else:
            rc_Qarea_min, rc_Qarea_max = 70*70, 90*90
        #
    elif irk==160:
        rc_d_ss = 156.
        rcFracOverlapOK = 0.25
        rc_maxDevMeanAreaQuads = 32
        rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf = 1.e-5, 1.e-5, 1.e-5
        rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.1, 0.1, 0.1
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig = 7.5e-2, 7.5e-2, 7.5e-2, 0.02
        if mode=='model':            
            rc_Qarea_min, rc_Qarea_max = 156*156, 164*164
        else:
            rc_Qarea_min, rc_Qarea_max = 140*140, 180*180
        #
    elif irk==320:
        rc_d_ss = 315.6
        rcFracOverlapOK = 0.25
        rc_maxDevMeanAreaQuads = 64
        rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf = 1.e-5, 1.e-5, 1.e-5
        rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.1, 0.1, 0.1
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig = 5.e-2, 5.e-2, 5.e-2, 0.01
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max =  312*312, 328*328
        else:
            rc_Qarea_min, rc_Qarea_max =  280*280, 360*360            
        #
    elif irk==640:
        rc_d_ss = 636.
        rcFracOverlapOK = 0.4
        rc_maxDevMeanAreaQuads = 128
        rc_div_min_pdf, rc_shr_min_pdf, rc_tot_min_pdf = 1.e-5, 1.e-5, 1.e-5
        rc_div_max_pdf, rc_shr_max_pdf, rc_tot_max_pdf =  0.1, 0.1, 0.1
        rc_div_max_fig, rc_shr_max_fig, rc_tot_max_fig, rc_df_fig = 2.5e-2, 2.5e-2, 2.5e-2, 0.01
        if mode=='model':
            rc_Qarea_min, rc_Qarea_max =  624*624, 656*656
        else:
            rc_Qarea_min, rc_Qarea_max =  480*480, 800*800
        #
    else:
        print('ERROR [updateConfig4Scale()]: scale "'+str(irk)+' km" is unknown!'); sys.exit(0)

    
    if mode in ['rgps','rgps_track']:
        from math import log2
        # scales = np.array([2**i * 10 for i in range(7)])
        # i      = np.array([ log2(rr/10) for rr in scales ])
        zdA_exp = 2.2
        zxp = int(log2(irk/10))
        zdA_tol = round( zdA_exp**zxp, 1 )
        rc_maxDevMeanAreaQuads = zdA_tol ; # this is lose 
        zrk = float(irk)
        if irk==10:
            zrk = 10.25
        ilw, iuw = zrk-zdA_tol , zrk+zdA_tol        
        rc_Qarea_min, rc_Qarea_max = ilw*ilw, iuw*iuw

    
    if ltalk: print(' *** [updateConfig4Scale](): upper and lower bound for scale "'+str(res_km)+' km":', np.sqrt([rc_Qarea_min, rc_Qarea_max]))
    if ltalk:
        print('              => rc_Qarea_min, Qarea_Nom, rc_Qarea_max =',rc_Qarea_min, res_km*res_km, rc_Qarea_max,'km^2')
        if mode in ['rgps','rgps_track']: print('                                   (+-',zdA_tol,'km)')

    #############################

    # Only for RGPS data:
    # *******************
    # When not at the nominal scale, we can adapt `rc_t_dev_cancel` to the scale we are dealing with:
    #    a quadrangle can involve points that are too distant from one another in terms of time
    #    => we disregard any quadrangles which standard deviation of the time of the 4 positions
    #       excess `t_dev_cancel` seconds !
    
    rc_t_dev_cancel = 60

    if mode in ['rgps','rgps_track','rgps_map']:

        #if res_km < 18:
        #    rc_t_dev_cancel = 60
        #else:
        #    rc_t_dev_cancel = 3*24*3600
            
        if 1==1:
            if   res_km>=15. and res_km<35.:
                rc_t_dev_cancel =  600
            elif   res_km>=35. and res_km<45.:
                rc_t_dev_cancel =  1200
                #
            elif   res_km>=70. and res_km<150.:
                rc_t_dev_cancel =  3600 ; #ok
            elif res_km>=150. and res_km<300.:
                rc_t_dev_cancel =  9*3600  ; #ok
            elif res_km>=300. and res_km<600.:
                rc_t_dev_cancel =  18*3600 ; #ok
            elif res_km>=600.:
                rc_t_dev_cancel = 24*3600 ; #ok...
                
    if rc_t_dev_cancel > 60:
        if ltalk: print('\n *** `rc_t_dev_cancel` updated to ',rc_t_dev_cancel/3600,'hours!')

    return 0
