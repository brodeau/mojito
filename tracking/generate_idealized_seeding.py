#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:

'''

from sys import argv, exit
from os import path, mkdir, makedirs, environ
import numpy as np
from re import split

import mojito   as mjt

idebug=0
iplot=1

# `seeding_type` is overriden if a NC file is specified (3rd argument) when calling this script...
#seeding_type='nemo_Tpoint' ; 
seeding_type='debug' ; # By default, ha

iHSS=15 ; # horizontal sub-sampling for NEMO points initialization...

ldo_coastal_clean = True
MinDistFromLand  = 100 ; # how far from the nearest coast should our buoys be? [km]
fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

if __name__ == '__main__':

    print('')
    if ldo_coastal_clean:
        cdata_dir = environ.get('DATA_DIR')
        if cdata_dir==None:
            print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
        fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'
    
    if not len(argv) in [2,3]:
        print('Usage: '+argv[0]+' <YYYY-MM-DD_hh:mm:ss> (<mesh_mask>,<iHSS> (for NEMO))')
        exit(0)

    cdate0 = argv[1]

    lNEMO = ( len(argv)==3 )
    if lNEMO:
        cvf = split(',',argv[2])
        cf_mm = cvf[0]
        seeding_type='nemo_Tpoint'
        iHSS = int(cvf[1])
        if iHSS<1 or iHSS>20:
            print('ERROR: chosen horizontal subsampling makes no sense iHSS=',iHSS)
            exit(0)

        
    

    if lNEMO:
        # Getting model grid metrics and friends:
        imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM = mjt.GetModelGrid( cf_mm )

        # A fake sea-ice concentration:
        #xIC[:,:] = id_uv.variables['siconc'][0,:,:] ; # We need ice conc. at t=0 so we can cancel buoys accordingly
        xIC = np.ones( np.shape(imaskt) )

    ############################
    # Initialization / Seeding #
    ############################


    if seeding_type=='nemo_Tpoint':
        XseedG = mjt.nemoSeed( imaskt, xlatT, xlonT, xIC, khss=iHSS, fmsk_rstrct=None )
        #
    elif seeding_type=='debug':
        if idebug in [0,1,2]: XseedG = mjt.debugSeeding()
        if idebug in [3]:     XseedG = mjt.debugSeeding1()
        #
    else:
         print(' ERROR: `seeding_type` =',seeding_type,' is unknown!')
         exit(0)
         
    print('\n * Shape of XseedG =',np.shape(XseedG))
    
    (nP,_) = np.shape(XseedG)




    zIDs = np.array([ i+1 for i in range(nP) ] , dtype=int)

    zTime = np.array( [ mjt.clock2epoch(cdate0) ], dtype='i4' )    
    print('\n * Requested initialization date =', mjt.epoch2clock(zTime[0]))
    
    print( zIDs  )
    print( zTime )

    if ldo_coastal_clean:
        mask = mjt.MaskCoastal( XseedG, rMinDistFromLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc, convArray='C' )
        print(' * Need to remove '+str(nP-np.sum(mask))+' points because too close to land! ('+str(MinDistFromLand)+'km)')
        nP = np.sum(mask) ; # new size once buoys too close to land removed...
        idxKeep    = np.where(mask==1)
        zIDs = zIDs[idxKeep]
        XseedG =  np.array( [ np.squeeze(XseedG[idxKeep,0]), np.squeeze(XseedG[idxKeep,1]) ] ).T
        del idxKeep
        
    XseedC = mjt.Geo2CartNPSkm1D( XseedG ) ; # same for seeded initial positions, XseedG->XseedC
    
    XseedG = np.reshape( XseedG, (1,nP,2) )
    XseedC = np.reshape( XseedC, (1,nP,2) )
    
    makedirs( './nc', exist_ok=True )
    foutnc = './nc/mojito_seeding_'+seeding_type+'_'+mjt.epoch2clock(zTime[0], precision='D')+'_HSS'+str(iHSS)+'.nc'
    
    print('\n *** Saving seeding file for date =',mjt.epoch2clock(zTime[0]))
    
    kk = mjt.ncSaveCloudBuoys( foutnc, zTime, zIDs, XseedC[:,:,0], XseedC[:,:,1], XseedG[:,:,0], XseedG[:,:,1],
                               corigin='idealized_seeding' )

    if iplot>0:

        makedirs( './figs', exist_ok=True )
        ffig = './figs/mojito_seeding_'+seeding_type+'_'+mjt.epoch2clock(zTime[0], precision='D')+'.png'

        
        mjt.ShowBuoysMap( 0, XseedG[0,:,1], XseedG[0,:,0], pvIDs=zIDs,
                          cfig=ffig, cnmfig=None, ms=5, ralpha=0.5, lShowDate=True,
                          zoom=1., title='Seeding initialization' )

    

