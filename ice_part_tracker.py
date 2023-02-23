#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:
   
 * `idateStop`: it must be possible to specify it via command line!
    => but if not specified + we are seeding from a netCDF RGPS data spawned by mojito
       ==> then get `idateStop` and `idateStart` from the seeding NC file !!!

 * So use a proper argument parsing stuff here
 
 * Now, when a `idateStart` exists, we have to locate the right time record to start from
   in the input SI3 UV file !!!!

 * The `Initialized_buoys_NANUK4_ICE-EVP00_1h_19970101_19970331_icemod.npz` file must bear the appropriate dates
   it's wrong now (it's that of the time span of the input SI3 UV file...)

'''

from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

from netCDF4 import Dataset

from climporn import dump_2d_field, epoch2clock, clock2epoch
import mojito   as mjt

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from math import atan2,pi

idebug=1
iplot=1

rdt = 3600. ; # time step (must be that of model output ice velocities used)

toDegrees = 180./pi

isubsamp_fig = 72 ; # frequency, in number of model records, we spawn a figure on the map (if idebug>2!!!)

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!

# `seeding_type` is overriden if a NC file is specified (3rd argument) when calling this script...
seeding_type='nemo_Tpoint' ; iHSS=15
#seeding_type='debug'

FillValue = -9999.


if __name__ == '__main__':

    if not len(argv) in [4]:
        print('Usage: '+argv[0]+' <ice_file_velocities_SI3> <mesh_mask> (<seeding_file.nc>)')
        exit(0)
        
    cf_uv = argv[1]
    cf_mm = argv[2]
    #cf_de = argv[3]

    idateStop = None
    
    lSeedFromNCfile = ( len(argv)==4 )
    if lSeedFromNCfile:
        seeding_type='mojitoNC'
        fNCseed = argv[3]
        mjt.chck4f(fNCseed)
        print('\n *** Will read initial seeding positions in first record of file:\n      => '+fNCseed+' !')
        ntr, zt = mjt.LoadNCtimeMJT( fNCseed, iverbose=idebug )
        idate0 = zt[0]     ; cdate0 = epoch2clock(idate0)
        idateN = zt[ntr-1] ; cdateN = epoch2clock(idateN)
        print('    => earliest and latest time position in the file: '+cdate0+' - '+cdateN)

        idate0 =  int( round( idate0/3600., 0 ) * 3600. ) ; cdate0 = epoch2clock(idate0)
        idateN =  int( round( idateN/3600., 0 ) * 3600. ) ; cdateN = epoch2clock(idateN)
        print('    ==> will actually use rounded to the hour! => '+cdate0+' - '+cdateN)

        idateStrt, cdateStrt = idate0, cdate0
        idateStop, cdateStop = idateN, cdateN

    #exit(0)
    #cdateStop = cf_de[0:4]+'-'+cf_de[4:6]+'-'+cf_de[6:8]+'_00:00:00'
    #idateStop = int(clock2epoch(cdateStop))
        
    cfdir = './figs/tracking'
    if iplot>0 and not path.exists(cfdir): mkdir(cfdir)
    if not path.exists('./npz'):           mkdir('./npz')


    # Infer name of NEMO CONFIG and experiment from SI3 file:
    vn = split('_',path.basename(cf_uv))
    nconf, nexpr = vn[0], split('-',vn[1])[1]
    print('    * NEMO config and experiment =', nconf, nexpr)
    
    # Getting model grid metrics and friends:
    imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM = mjt.GetModelGrid( cf_mm )

    # Allocating arrays for model data:
    (Nj,Ni) = np.shape( imaskt )    
    xUu = np.zeros((Nj,Ni))
    xVv = np.zeros((Nj,Ni))
    xIC = np.zeros((Nj,Ni)) ; # Sea-ice concentration

    # We need a name for the intermediate backup file:
    if seeding_type=='nemo_Tpoint':
        cf_npz_intrmdt = './npz/Initialized_buoys_'+seeding_type+'_HSS'+str(iHSS)+'.npz'
    elif seeding_type=='debug':
        cf_npz_intrmdt = './npz/Initialized_buoys_'+seeding_type+'.npz'
    elif seeding_type=='mojitoNC':
        cf_npz_intrmdt = './npz/Initialized_buoys_'+str.replace( path.basename(cf_uv), '.nc4', '.npz' )
    else:
        print('ERROR: "'+seeding_type+'" is an unknown seeding option!')
        exit(0)

    
    # Open Input SI3 data file
    id_uv = Dataset(cf_uv)
    Nt = id_uv.dimensions['time_counter'].size
    if id_uv.variables['time_counter'].units != ctunits_expected:
        print('ERROR: wrong units for time calendar in file:',cf_uv)
        exit(0)
    print('\n\n *** '+str(Nt)+' records in input SI3 file!')
    ztime_model = np.array( id_uv.variables['time_counter'][:] , dtype='i4' )


    if not idateStop:
        print('FIXME: we need a date stop `idateStop` !'); exit(0)
    
    if idateStop<=np.min(ztime_model) or idateStop>np.max(ztime_model):
        print('ERROR: end time is outside of time calendar of the model file!')
        exit(0)
    kstop = np.argmin(np.abs(ztime_model[:]-idateStop)) + 1
    print('   => will stop at record',kstop,' =>',epoch2clock(ztime_model[kstop]))
    Nt = kstop

    
    ############################
    # Initialization / Seeding #
    ############################


    if path.exists(cf_npz_intrmdt):

        print('\n *** We found file '+cf_npz_intrmdt+' here! So using it and skeeping first stage!')
        with np.load(cf_npz_intrmdt) as data:
            nP    = data['nP']
            xPosG0 = data['xPosG0']
            xPosC0 = data['xPosC0']
            IDs   = data['IDs']
            vJIt  = data['vJIt']
            VRTCS = data['VRTCS']

    else:

        # Going through whole initialization / seeding process
        # ----------------------------------------------------
                
        xIC[:,:] = id_uv.variables['siconc'][0,:,:] ; # We need ice conc. at t=0 so we can cancel buoys accordingly
    
        if seeding_type=='nemo_Tpoint':
            XseedG = mjt.nemoSeed( imaskt, xlatT, xlonT, xIC, khss=iHSS, fmsk_rstrct=None )
            #
        elif seeding_type=='debug':
            if idebug in [0,1,2]: XseedG = mjt.debugSeeding()
            if idebug in [3]:     XseedG = mjt.debugSeeding1()
            #
        elif seeding_type=='mojitoNC':
            zt, zIDs, XseedG, XseedC = mjt.LoadNCdataMJT( fNCseed, krec=0, iverbose=idebug )
            print('     => data is read at date =',epoch2clock(zt))
        #
        print('\n shape of XseedG =',np.shape(XseedG))
    
        (nP,_) = np.shape(XseedG)

        # We want an ID for each seeded buoy:
        IDs = np.array( range(nP), dtype=int) + 1 ; # Default! No ID=0 !!!

        
        if seeding_type == 'mojitoNC':
            IDs[:] = zIDs[:]
            del zIDs
        else:
            XseedC = mjt.Geo2CartNPSkm1D( XseedG ) ; # same for seeded initial positions, XseedG->XseedC
            
        # Find the location of each seeded buoy onto the model grid:
        nP, xPosG0, xPosC0, IDs, vJIt, VRTCS = mjt.SeedInit( IDs, XseedG, XseedC, xlatT, xlonT, xYf, xXf,
                                                             xResKM, imaskt, xIceConc=xIC, iverbose=idebug )    
        del XseedG, XseedC
        
        # This first stage is fairly costly so saving the info:
        print('\n *** Saving intermediate data into '+cf_npz_intrmdt+'!')
        np.savez_compressed( cf_npz_intrmdt, nP=nP, xPosG0=xPosG0, xPosC0=xPosC0, IDs=IDs, vJIt=vJIt, VRTCS=VRTCS )

    del xResKM

        
    # Allocation for nP buoys:
    iAlive = np.zeros(      nP , dtype='i1') + 1 ; # tells if a buoy is alive (1) or zombie (0) (discontinued)
    vTime  = np.zeros( Nt+1, dtype=int ) ; # UNIX epoch time associated to position below
    xmask  = np.zeros((Nt+1,nP,2), dtype='i1')
    xPosC  = np.zeros((Nt+1,nP,2)) + FillValue  ; # x-position of buoy along the Nt records [km]
    xPosG  = np.zeros((Nt+1,nP,2)) + FillValue
    vCELLs   = np.zeros(nP, dtype=Polygon) ; # stores for each buoy the polygon object associated to the current mesh/cell
    lStillIn = np.zeros(nP, dtype=bool) ; # tells if a buoy is still within expected mesh/cell..

    # Initial values for some arrays:
    xPosC[0,:,:] = xPosC0
    xPosG[0,:,:] = xPosG0    
    xmask[0,:,:] = 1
    
    del xPosC0, xPosG0

    if iplot>0 and idebug>0:
        mjt.ShowBuoysMap( 0, xPosG[0,:,1], xPosG[0,:,0], pvIDs=IDs, cfig=cfdir+'/INIT_Pos_buoys_'+'%4.4i'%(0)+'.png',
                          cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1., title='IceTracker: Init Seeding' )

    
    ######################################
    # Loop along model data time records #
    ######################################

    for jt in range(Nt):
        rtmod = id_uv.variables['time_counter'][jt] ; # time of model data (center of the average period which should = rdt)
        itime = int(rtmod - rdt/2.) ; # velocitie is average under the whole rdt, at the center!
        ctime = epoch2clock(itime)
        print('\n *** Reading record #'+str(jt)+'/'+str(Nt-1)+' in SI3 file ==> date =',
              ctime,'(model:'+epoch2clock(int(rtmod))+')')
        vTime[jt] = itime
        
        xIC[:,:] = id_uv.variables['siconc'][jt,:,:]
        xUu[:,:] = id_uv.variables['u_ice'][jt,:,:]
        xVv[:,:] = id_uv.variables['v_ice'][jt,:,:]

        print('   *   current number of buoys alive = '+str(iAlive.sum()))
        
        for jP in range(nP):

            if iAlive[jP]==1:
            
                [ ry  , rx   ] = xPosC[jt,jP,:] ; # km !
                [ rlat, rlon ] = xPosG[jt,jP,:] ; # degrees !
    
                [jnT,inT] = vJIt[jP,:]; # indices for T-point at center of current cell
                
                if idebug>0:
                    print('\n    * BUOY ID:'+str(IDs[jP])+' => jt='+str(jt)+': ry, rx =', ry, rx,'km'+': rlat, rlon =', rlat, rlon )
    
                ######################### N E W   M E S H   R E L O C A T I O N #################################
                if not lStillIn[jP]:
                    if idebug>0: print('      +++ RELOCATION NEEDED for buoy with ID:'+str(IDs[jP])+' +++')
    
                    [ [ jbl, jbr, jur, jul ], [ ibl, ibr, iur, iul ] ] = VRTCS[jP,:,:]

                    if idebug>0:
                        print('     ==> 4 corner points of our mesh (anti-clockwise, starting from BLC) =',
                              [ [jbl,ibl],[jbr,ibr],[jur,iur],[jul,iul] ])
    
                    # The current mesh/cell as a shapely polygon object:
                    vCELLs[jP] = Polygon( [ (xYf[jbl,ibl],xXf[jbl,ibl]) , (xYf[jbr,ibr],xXf[jbr,ibr]) ,
                                            (xYf[jur,iur],xXf[jur,iur]) , (xYf[jul,iul],xXf[jul,iul]) ] )
                    
                    if idebug>0:
                        # Control if we are dealing with the proper cell/mesh
                        #   * buoy location ?
                        if not mjt.IsInsideCell(ry, rx, vCELLs[jP]):
                            print('\nPROBLEM: buoy location is not inside the expected mesh!!!')
                            print([(xYf[jbl,ibl],xXf[jbl,ibl]), (xYf[jbr,ibr],xXf[jbr,ibr]), (xYf[jur,iur],xXf[jur,iur]),
                                   (xYf[jul,iul],xXf[jul,iul])])
                            exit(0)
                        #
                        #    * T-point @ center of mesh ?
                        if not mjt.IsInsideCell( xYt[jnT,inT], xXt[jnT,inT], vCELLs[jP]):
                            print('\nPROBLEM: T-point is not inside the expected mesh!!!')
                            print([(xYf[jbl,ibl],xXf[jbl,ibl]), (xYf[jbr,ibr],xXf[jbr,ibr]), (xYf[jur,iur],xXf[jur,iur]),
                                   (xYf[jul,iul],xXf[jul,iul])])
                            exit(0)
    
                        print('     +++ control of location of point inside selected cell successfuly passed! :D')
                        #if jt>0: exit(0); #fixme rm!
                ### if not lStillIn[jP]
                ###########################################################################################################################
    
                if idebug>2:
                    # We can have a look in the mesh:
                    cnames    = np.array([ 'P'+str(i+1)+': '+str(VRTCS[jP,0,i])+','+str(VRTCS[jP,1,i]) for i in range(4) ], dtype='U32')
                    #mjt.PlotMesh( (rlat,rlon), xlatF, xlonF, VRTCS[jP,:,:].T, vnames=cnames,
                    #              fig_name='mesh_lon-lat_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                    #              pcoor_extra=(xlatT[jnT,inT],xlonT[jnT,inT]), label_extra='T-point' )
                    mjt.PlotMesh( ( ry , rx ),  xYf ,  xXf,  VRTCS[jP,:,:].T, vnames=cnames,
                                  fig_name='mesh_X-Y_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                                  pcoor_extra=(xYt[jnT,inT],xXt[jnT,inT]), label_extra='T-point' )
    

                # j,i indices of the cell we are dealing with = that of the F-point aka the upper-right point !!!
                [jT,iT] = vJIt[jP,:]
                
                # ASSUMING THAT THE ENTIRE CELL IS MOVING AT THE SAME VELOCITY: THAT OF U-POINT OF CELL
                # zU, zV = xUu[jT,iT], xVv[jT,iT] ; # because the F-point is the upper-right corner
                zU = 0.5*(xUu[jT,iT]+xUu[jT,iT-1])
                zV = 0.5*(xVv[jT,iT]+xVv[jT-1,iT])                
                if idebug>0:
                    print('    =>> read velocity at ji,jj=',iT,jT)
                    print('    * ice velocity of the mesh: u,v =',zU, zV, 'm/s')
    
                # Displacement during the upcomming time step:
                dx = zU*rdt
                dy = zV*rdt
                if idebug>0: print('      ==> displacement during `dt`: dx,dy =',dx,dy, 'm')

                # => position [km] of buoy after this time step will be:
                rx_nxt = rx + dx/1000. ; # [km]
                ry_nxt = ry + dy/1000. ; # [km]
                xPosC[jt+1,jP,:] = [ ry_nxt, rx_nxt ]
                xmask[jt+1,jP,:] = [    1  ,    1   ]
                
                # Is it still inside our mesh:
                lSI = mjt.IsInsideCell(ry_nxt, rx_nxt, vCELLs[jP])
                lStillIn[jP] = lSI
                if idebug>0: print('      ==> Still inside the same mesh???',lSI)
    
                if not lSI:
                    # => point is exiting the current cell!
                    
                    # Tells which of the 4 cell walls the point has crossed:
                    icross = mjt.CrossedEdge( [ry,rx], [ry_nxt,rx_nxt], VRTCS[jP,:,:], xYf, xXf, iverbose=idebug )
                    
                    # Tells in which adjacent cell the point has moved:
                    inhc = mjt.NewHostCell( icross, [ry,rx], [ry_nxt,rx_nxt], VRTCS[jP,:,:], xYf, xXf,  iverbose=idebug )
                    
                    # Update the mesh indices according to the new host cell:
                    VRTCS[jP,:,:],vJIt[jP,:] = mjt.UpdtInd4NewCell( inhc, VRTCS[jP,:,:], vJIt[jP,:] )
                    # LOLO DEBUG: test if `UpdtInd4NewCell` does a good job:
                    #[ [ jbl, jbr, jur, jul ], [ ibl, ibr, iur, iul ] ] = VRTCS[jP,:,:]
                    #zzf = Polygon( [ (xYf[jbl,ibl],xXf[jbl,ibl]) , (xYf[jbr,ibr],xXf[jbr,ibr]) ,
                    #                        (xYf[jur,iur],xXf[jur,iur]) , (xYf[jul,iul],xXf[jul,iul]) ] )
                    #if not mjt.IsInsideCell( ry_nxt, rx_nxt, zzf ):
                    #    print('FUCK UP!!!')
                    #    exit(0)
                    #else:
                    #    print('INSIDE cell :D')
                                        
                    # Based on updated new indices, some buoys might get killed:
                    icncl = mjt.Survive( IDs[jP], vJIt[jP,:], imaskt, pIceC=xIC,  iverbose=idebug )
                    if icncl>0: iAlive[jP]=0
                    
                ### if not lSI

            ### if iAlive[jP]==1
                
        ### for jP in range(nP)

        # Updating in terms of lon,lat for all the buoys at once:
        xPosG[jt+1,:,:] = mjt.CartNPSkm2Geo1D( xPosC[jt+1,:,:] )


        if iplot>0:
            # Show on the map of the Arctic:
            if jt%isubsamp_fig == 0:
                zLon = np.ma.masked_where( xmask[jt,:,1]==0, xPosG[jt,:,1] )
                zLat = np.ma.masked_where( xmask[jt,:,0]==0, xPosG[jt,:,0] )
                mjt.ShowBuoysMap( itime, zLon, zLat, pvIDs=IDs,
                                  cfig=cfdir+'/Pos_buoys_'+nexpr+'_'+'%4.4i'%(jt)+'_'+ctime+'.png',
                                  cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1.,
                                  title='IceTracker + SI3 '+nexpr+' u,v fields' )
                del zLon, zLat
        print('\n')                
    ### for jt in range(Nt-1)

    id_uv.close()

    vTime[Nt] = vTime[Nt-1] + int(rdt)


    # Masking arrays:
    xPosG = np.ma.masked_where( xmask==0, xPosG )
    xPosC = np.ma.masked_where( xmask==0, xPosC )
    
    # ==> time to save itime, xPosXX, xPosYY, xPosLo, xPosLa into a netCDF file !
    cdt1 = split('_', epoch2clock(vTime[0] ))[0]
    cdt2 = split('_', epoch2clock(vTime[Nt]))[0]
    corgn = 'NEMO-SI3_'+nconf+'_'+nexpr
    foutnc = corgn+'_tracking_'+cdt1+'_'+cdt2+'.nc'
    
    kk = mjt.ncSaveCloudBuoys( foutnc, vTime, IDs, xPosC[:,:,0], xPosC[:,:,1], xPosG[:,:,0], xPosG[:,:,1],
                               mask=xmask[:,:,0], tunits=ctunits_expected, fillVal=FillValue, corigin=corgn )
            
