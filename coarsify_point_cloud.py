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
from os import path, environ, mkdir, makedirs
import numpy as np
from re import split
from scipy.spatial import Delaunay

import mojito   as mjt
from mojito import config as cfg

idebug = 0
iplot  = 1 ; # Create figures to see what we are doing...

rzoom_fig = 2

l_varying_dist2coast = False
l_intermediate_randomization_big_scales = True

if __name__ == '__main__':

    #cdata_dir = environ.get('DATA_DIR')

    
    rd_ss = None
    rd_tc = None
    
    if not len(argv) in [4,5,6]:
        print('Usage: '+argv[0]+' <mode:rgps/model> <file_mojito.nc> <res_km> (<rd_ss_km>) (<min_dist_coast>)')
        exit(0)

    # Idea remove all the buoys closer to land than `rd_ss_km` km !!!
    # => this way avoids big quadrangles attached to coastal regions
    # => better for randomizing...    

    quality_mode  = argv[1]
    cf_nc_in = argv[2]
    creskm = argv[3]
    reskm = float(creskm)

    ik = cfg.controlModeName( path.basename(__file__), quality_mode )

    kk = cfg.initialize( mode=quality_mode )
        
    ldss = ( len(argv)==5 or len(argv)==6 ); # a `rd_ss` is provided!
    ldtc =   len(argv)==6 ; # a minimum distance to coast is provided

    if ldss:
        rd_ss = float(argv[4])
    else:
        kl = cfg.updateConfig4Scale( int(reskm), mode=quality_mode )
        rd_ss = cfg.rc_d_ss
        
    if ldtc:        
        rd_tc = float(argv[5])
        
    mjt.chck4f(cf_nc_in)

    if not path.exists('./nc'): mkdir('./nc')
    if iplot>0:
        fdir = './figs/coarsify/'+creskm+'km'
        makedirs( fdir, exist_ok=True )

    #########################################################################################################



    cfbn = str.replace( path.basename(cf_nc_in), '.nc', '')
    if ldss:
        cf_nc_out = './nc/'+cfbn+'_'+str(int(rd_ss))+'-'+creskm+'km.nc'
    else:
        cf_nc_out = './nc/'+cfbn+'_'+creskm+'km.nc'
    cfbase = str.replace( cfbn, 'SELECTION_', '')
    print('     => will generate: '+cf_nc_out,'\n')


    if not path.exists(cf_nc_out):
            
        print(' *** Will coarsify for a spatial scale of '+creskm+' km')
    
        
        # Getting dimensions and some info:
        Nrec, nBmax, corigin, lTimePos = mjt.GetDimNCdataMJT( cf_nc_in )
    
        if lTimePos: print(' *** There is the "time_pos" data in input netCDF file! => gonna use it!')
    
        vdate = np.zeros( Nrec,  dtype=int )
        vIDs  = np.zeros( nBmax, dtype=int )
        zGC   = np.zeros( (nBmax,2,Nrec) )
        zXY   = np.zeros( (nBmax,2,Nrec) )
        zmsk  = np.zeros( (nBmax,Nrec), dtype='i1' )
        ztim  = np.zeros( (nBmax,Nrec), dtype=int )        

        if lTimePos:        
            vdate[:], vIDs[:], zGC[:,:,:], zXY[:,:,:], zmsk[:,:], ztim[:,:] = mjt.LoadNCdataMJT( cf_nc_in, lmask=True, lGetTimePos=True,  convention='F' )
        else:
            vdate[:], vIDs[:], zGC[:,:,:], zXY[:,:,:], zmsk[:,:]            = mjt.LoadNCdataMJT( cf_nc_in, lmask=True, lGetTimePos=False, convention='F' )
    
        # Need some calendar info:
        NbDays = int( (vdate[1] - vdate[0]) )/ cfg.rc_day2sec
        cdt1 = mjt.epoch2clock(vdate[0] )
        cdt2 = mjt.epoch2clock(vdate[-1])
    
        print('    *  start and End dates => '+cdt1+' -- '+cdt2,' | number of buoys =>',np.sum(zmsk[:,0]), np.sum(zmsk[:,1]))
        print('        ==> nb of days =', NbDays)
    
        mask = np.zeros( nBmax      , dtype='i1') + 1  ; # Mask to for "deleted" points (to cancel)    

        zPnm = np.array( [ str(i) for i in vIDs ], dtype='U32' ) ; # Name for each point, based on 1st record...
    
        if not lTimePos:
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
    

        
        if l_varying_dist2coast and ldss:

            rDmin = 100.
            
            if   ldtc:
                rDmin = rd_tc
            elif corigin=='RGPS' and reskm>300:
                from random import random
                zrnd = 2*random() - 1 ; # Between -1 and 1
                rDmin =  150. - 50.*zrnd
                print('    LOLO: rDmin =',rDmin)
            #
            #
            if rDmin>101.:
                print( ' *** COARSIFICATION: reskm=',reskm,'=> rd_ss=',rd_ss,'=> MaskCoastal() with `rDmin` =',rDmin,'km !!!' )
                
                for jr in range(Nrec):
                    print('    * Record #'+str(jr)+':')
                    mask[:] = mjt.MaskCoastal( zGC[:,:,jr], mask=mask[:], rMinDistLand=rDmin, fNCdist2coast=cfg.fdist2coast_nc, convArray='F' )
                # How many points left after elimination of buoys that get too close to land (at any record):
                NbP  = np.sum(mask)
                NbRM = nBmax-NbP
                print('\n *** '+str(NbP)+' / '+str(nBmax)+' points survived the `dist2coast` cleanup => ', str(NbRM)+' points deleted!')
                if NbRM>0:
                    zPnm, vIDs, zGC, zXY, ztim = mjt.ShrinkArrays( mask, zPnm, vIDs, zGC, zXY, ztim, recAxis=2 )
                if NbP<4:
                    print('\n *** Exiting because no enough points alive left!')
                    exit(0)
            
        ### if l_varying_dist2coast and ldss




                    
        # We must call the coarsening function only for first record !
        # => following records are just the same coarsened buoys...
        
    
        jr = 0
    
        print('\n\n *** COARSIFICATION => record jr='+str(jr))
    
        cdats  = mjt.epoch2clock(vdate[jr])
        cdate  = str.replace( mjt.epoch2clock(vdate[jr], precision='D'), '-', '')
    
        print('    * which is original record '+str(jr)+' => date =',cdats)


        if ldss and corigin=='RGPS' and reskm>300 and l_intermediate_randomization_big_scales:
            #
            # We do it in a 2-stage fashion introcuding noise
            # A/ to something betwee 8 and 30 km:
            zrnd = 2*random() - 1 ; # Between -1 and +1
            rd_int = 0.5*( (8. + 30.) - (30. - 8.)*zrnd )
            print(' ### Intermediate coarsening radius =',rd_int,'km')
            NbPss0, zXYss0, idxKeep0 = mjt.SubSampCloud( rd_int, zXY[:,:,jr] )
            print(' ### Now, final/post-intermediate coarsening radius =',rd_ss,'km')
            NbPss, zXYss, idxKeep2 = mjt.SubSampCloud( rd_ss, zXYss0[:,:] )
            idxKeep = idxKeep0[idxKeep2]
            del NbPss0, zXYss0, idxKeep0, idxKeep2
            #
        else:
            print('\n *** Applying spatial sub-sampling with radius: '+str(round(rd_ss,2))+'km for record jr=',jr)
            NbPss, zXYss, idxKeep = mjt.SubSampCloud( rd_ss, zXY[:,:,jr] )
    
        zYkm = np.zeros( (Nrec,NbPss) )
        zXkm = np.zeros( (Nrec,NbPss) )
        zlat = np.zeros( (Nrec,NbPss) )
        zlon = np.zeros( (Nrec,NbPss) )
    
        zYkm[:,:] = zXY[idxKeep,1,:].T
        zlat[:,:] = zGC[idxKeep,1,:].T    
        zXkm[:,:] = zXY[idxKeep,0,:].T
        zlon[:,:] = zGC[idxKeep,0,:].T
        
        print('   * saving '+cf_nc_out)
    
        kk = mjt.ncSaveCloudBuoys( cf_nc_out, vdate, vIDs[idxKeep], zYkm, zXkm, zlat, zlon, mask=zmsk[idxKeep,:].T,
                                   xtime=ztim[idxKeep,:].T, fillVal=mjt.FillValue, corigin=corigin )
    
        print()
    
    
        # Some debug plots to see the job done:
    
        if iplot>0:
            print('\n *** Some plots...')
            
            for jr in range(Nrec):
            
                cfb = cfbase+'_rec%3.3i'%(jr)+'_rdss'+str(int(rd_ss))
                cfc = '_SS'+creskm+'km'
                if ldtc:
                    cfc += '_DCmin'+str(int(rd_tc))
    
                if iplot>1:
                    # A: on cartesian grid
                    #######################      
                    if jr==0:
                        # We need to find a descent X and Y range for the figures:
                        vrngX = mjt.roundAxisRange( zXY[:,0,0], rndKM=50. )
                        vrngY = mjt.roundAxisRange( zXY[:,1,0], rndKM=50. )
                    
                    if idebug>0:
                        kk = mjt.ShowTQMesh( zXY[:,0,jr], zXY[:,1,jr], cfig=fdir+'/00_'+cfb+'_Original.png',
                                             lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY ) ; #ppntIDs=vIDs[:], 
                    
                    # After subsampling
                    kk = mjt.ShowTQMesh( zXkm[jr,:], zYkm[jr,:], cfig=fdir+'/00_'+cfb+cfc+'.png',
                                         lGeoCoor=False, zoom=rzoom_fig, rangeX=vrngX, rangeY=vrngY ) ; #ppntIDs=vIDs[idxKeep], 
    
                
                # B: same, but on projection with geographic coordinates
                ########################################################
                if idebug>0:
                    kk = mjt.ShowBuoysMap( vdate[jr], zGC[jr,:,0], zGC[jr,:,1], pvIDs=vIDs[:], cfig=fdir+'/00_PROJ_'+cfb+'_Original.png',
                                           cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )
    
                kk = mjt.ShowBuoysMap( vdate[jr], zlon[jr,:], zlat[jr,:], pvIDs=vIDs[:], cfig=fdir+'/00_PROJ_'+cfb+cfc+'.png',
                                       cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )


    
    else:
        print('\n *** File "'+cf_nc_out+'" already exists!!!\n')
        
