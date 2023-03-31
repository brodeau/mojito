#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#  INPUT DATA: a `nc` file created with `ncSaveCloudBuoys()` of `mojito/ncio.py` !
#
#  OUTPUT: similar netcdf file but with less buoys (due to coarsening)
#
#    L. Brodeau, 2023
#
##################################################################

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np
from re import split
from scipy.spatial import Delaunay

import mojito   as mjt

idebug = 0
iplot  = 2 ; # Create figures to see what we are doing...

rzoom_fig = 2

if __name__ == '__main__':

    #cdata_dir = environ.get('DATA_DIR')

    if len(argv) != 3:
        print('Usage: '+argv[0]+' <file_mojito.nc> <res_km>')
        exit(0)

    cf_nc_in = argv[1]
    creskm = argv[2]
    reskm = float(creskm)

    #lNeed4SubSamp = ( reskm > 12. )
    
    mjt.chck4f(cf_nc_in)

    if not path.exists('./nc'): mkdir('./nc')
    if iplot>0:
        cfdir = './figs/coarsify'
        if not path.exists(cfdir): mkdir(cfdir)

    #########################################################################################################


    print(' *** Will coarsify for a spatial scale of '+creskm+' km')

    cfbn = str.replace( path.basename(cf_nc_in), '.nc', '')
    cf_nc_out = './nc/'+cfbn+'_'+creskm+'km.nc'
    cfbase = str.replace( cfbn, 'SELECTION_', '')
    print('     => will generate: '+cf_nc_out,'\n')

    
    # Getting dimensions and some info:
    Nrec, nBmax, corigin, lTimePos = mjt.GetDimNCdataMJT( cf_nc_in )

    if lTimePos: print(' *** There is the "time_pos" data in input netCDF file! => gonna use it!')



    #
    vdate = np.zeros( Nrec,  dtype=int )
    vIDs  = np.zeros( nBmax, dtype=int )
    xPosG = np.zeros( (Nrec,nBmax,2) )
    xPosC = np.zeros( (Nrec,nBmax,2) )
    pmsk  = np.zeros( (Nrec,nBmax), dtype='i1' )
    if lTimePos:
        timePos = np.zeros( (Nrec,nBmax), dtype=int )        
    #
    #for jr in range(Nrec):
    if lTimePos:        
        vdate[:], vIDs[:], xPosG[:,:,:], xPosC[:,:,:], pmsk[:,:], timePos[:,:] = mjt.LoadNCdataMJT( cf_nc_in, lmask=True, lGetTimePos=True  )
        #jr = 1
        #vdate[jr], vIDs[:], xPosG[jr,:,:], xPosC[jr,:,:], pmsk[jr,:], timePos[jr,:] = mjt.LoadNCdataMJT( cf_nc_in, krec=jr, lmask=True, lGetTimePos=True  )
    else:
        vdate[:], vIDs[:], xPosG[:,:,:], xPosC[:,:,:], pmsk[:,:]                = mjt.LoadNCdataMJT( cf_nc_in, lmask=True, lGetTimePos=False )
    #
    #print('LOLO: shape of vdate:', np.shape(vdate))
    #print('LOLO: shape of vIDs:', np.shape(vIDs))
    #print('LOLO: shape of xPosG:', np.shape(xPosG))
    #print('LOLO: shape of xPosC:', np.shape(xPosC))
    #print('LOLO: shape of pmsk:', np.shape(pmsk))
    #print('LOLO: shape of timePos:', np.shape(timePos))
    
    #print( ' * jr = ',jr, ', mean date =',mjt.epoch2clock(vdate[jr]))    
    #if jr==0:
    #    vIDs[:] = zIDs[:]
    #else:
    #    if np.sum(zIDs[:]-vIDs[:])!=0:
    #    print('ERROR: ID fuck up in input file!') ; exit(0)

    
    # Need some calendar info:
    NbDays = int( (vdate[1] - vdate[0]) / (3600.*24.) )
    cdt1 = mjt.epoch2clock(vdate[0] )
    cdt2 = mjt.epoch2clock(vdate[-1])

    print('    *  start and End dates => '+cdt1+' -- '+cdt2,' | number of buoys =>',np.sum(pmsk[0,:]), np.sum(pmsk[1,:]))
    print('        ==> nb of days =', NbDays)

    
    # STUPID: #fixme
    zXY   = np.zeros( (Nrec,nBmax,2) )
    zXY[:,:,0] = xPosC[:,:,1]
    zXY[:,:,1] = xPosC[:,:,0]
    zGC   = np.zeros( (Nrec,nBmax,2) )
    zGC[:,:,0] = xPosG[:,:,1]
    zGC[:,:,1] = xPosG[:,:,0]

    mask = np.zeros( nBmax      , dtype='i1') + 1  ; # Mask to for "deleted" points (to cancel)    
    ztim = np.zeros((Nrec,nBmax), dtype=int )
    zPnm = np.array( [ str(i) for i in vIDs ], dtype='U32' ) ; # Name for each point, based on 1st record...

    if lTimePos:
        ztim[:,:] = timePos[:,:]
    else:
        for jp in range(nBmax): ztim[jp,:] = vdate[:]


        
    NbP = nBmax



    if idebug>0:
        for jr in range(Nrec):
            print('\n  DEBUG *** Record jr='+str(jr)+':')
            for jc in range(NbP):
                print( ' * #'+str(jc)+' => Name: "'+zPnm[jc]+'": ID=',vIDs[jc],', X=',zXY[jc,0,jr],', Y=',zXY[jc,1,jr],
                       ', lon=',zGC[jc,0,jr],', lat=',zGC[jc,1,jr], ', time=',mjt.epoch2clock(ztim[jc,jr]) )
                if str(vIDs[jc])!=zPnm[jc]:
                    print(' Fuck Up!!!! => vIDs[jc], zPnm[jc] =',vIDs[jc], zPnm[jc] ); exit(0)
        print('')


    cdate0  = str.replace( mjt.epoch2clock(vdate[0], precision='D'), '-', '')





    # We must call the coarsening function only for first record !
    # => following records are just the same coarsened buoys...
    

    jr = 0

    print('\n\n *** COARSIFICATION => record jr='+str(jr))

    cdats  = mjt.epoch2clock(vdate[jr])
    cdate  = str.replace( mjt.epoch2clock(vdate[jr], precision='D'), '-', '')

    print('    * which is original record '+str(jr)+' => date =',cdats)


    # SUB-SAMPLING
    #  From experience we have to pick a scale significantly smaller that that of the desired
    #  quadrangles in order to obtain quadrangle of the correct size!

    if reskm==20.:
        rd_ss = 15.1
    elif reskm==40.:
        rd_ss = 35.
    elif reskm==80.:
        rd_ss = 73.
    elif reskm==160.:
        rd_ss = 145.
    elif reskm==320.:
        rd_ss = 295.
    elif reskm==640.:
        rd_ss = 620.
    else:
        print('ERROR: dont know what `rd_ss` to use for resolution ='+creskm+'km')
        exit(0)
        

    print('\n *** Applying spatial sub-sampling with radius: '+str(round(rd_ss,2))+'km for record jr=',jr)
    #print('LOLO: shape(zXY) =',np.shape(zXY))
    NbPss, zXYss, idxKeep = mjt.SubSampCloud( rd_ss, zXY[jr,:,:] )

    # vIDs[idxKeep], ztim[idxKeep,jr], zPnm[idxKeep]


    zYkm = np.zeros( (Nrec,NbPss) )
    zXkm = np.zeros( (Nrec,NbPss) )
    zlat = np.zeros( (Nrec,NbPss) )
    zlon = np.zeros( (Nrec,NbPss) )

    zYkm[:,:] = zXY[:,idxKeep,1]
    zlat[:,:] = zGC[:,idxKeep,1]    
    zXkm[:,:] = zXY[:,idxKeep,0]
    zlon[:,:] = zGC[:,idxKeep,0]

       
    print('   * saving '+cf_nc_out)

    kk = mjt.ncSaveCloudBuoys( cf_nc_out, vdate, vIDs[idxKeep], zYkm, zXkm, zlat, zlon, mask=pmsk[:,idxKeep],
                               xtime=ztim[:,idxKeep], fillVal=mjt.FillValue, corigin=corigin )

    print()


    # Some debug plots to see the job done:

    if iplot>0:
        print('\n *** Some plots...')
        
        for jr in range(Nrec):
        
            cfb = cfbase+'_rec%3.3i'%(jr)
            cfc = '_SS'+creskm+'km'

            # A: on cartesian grid
            #######################
            
            if jr==0:
                # We need to find a descent X and Y range for the figures:
                vrngX = mjt.roundAxisRange( zXY[0,:,0], rndKM=50. )
                vrngY = mjt.roundAxisRange( zXY[0,:,1], rndKM=50. )
                
            # Shows the cloud of buoys (with buoys' IDs) on the Cartesian plane (km)
            # Before subsampling
            kk = mjt.ShowTQMesh( zXY[jr,:,0], zXY[jr,:,1], cfig=cfdir+'/00_'+cfb+'_Original.png',
                                 lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY ) ; #ppntIDs=vIDs[:], 
            
            # After subsampling
            kk = mjt.ShowTQMesh( zXkm[jr,:], zYkm[jr,:], cfig=cfdir+'/00_'+cfb+cfc+'.png',
                                 lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY ) ; #ppntIDs=vIDs[idxKeep], 

            
            # B: same, but on projection with geographic coordinates
            ########################################################
            
            kk = mjt.ShowBuoysMap( vdate[jr], zGC[jr,:,0], zGC[jr,:,1], pvIDs=vIDs[:], cfig=cfdir+'/00_PROJ_'+cfb+'_Original.png',
                                   cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )

            kk = mjt.ShowBuoysMap( vdate[jr], zlon[jr,:], zlat[jr,:], pvIDs=vIDs[:], cfig=cfdir+'/00_PROJ_'+cfb+cfc+'.png',
                                   cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )





    
    #exit(0)
                
    
