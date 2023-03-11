#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:

   *

'''

from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

from netCDF4 import Dataset

import mojito   as mjt

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from math import atan2,pi

idebug=0
iplot=1

rdt = 3600 ; # time step (must be that of model output ice velocities used)

toDegrees = 180./pi

isubsamp_fig = 72 ; # frequency, in number of model records, we spawn a figure on the map (if idebug>2!!!)

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!

FillValue = -9999.

iUVstrategy = 1 ; #  What U,V should we use inside a given T-cell of the model?
#                 #  * 0 => use the same MEAN velocity in the whole cell => U = 0.5*(U[j,i-1] + U[j,i]), V = 0.5*(V[j-1,i] + U[j,i])
#                 #  * 1 => use the same NEAREST velocity in the whole cell => U = U[@ nearest U-point], V = V[@ nearest V-point]




if __name__ == '__main__':


    print('')
    print('##########################################################')
    print('#            MOJITO ICE PARTICULES TRACKER               #')
    print('##########################################################\n')

    if not len(argv) in [4,5]:
        print('Usage: '+argv[0]+' <ice_file_velocities_SI3> <mesh_mask> <seeding_file.nc> (<#record_to_seed_with>)')
        exit(0)

    cf_uv = argv[1]
    cf_mm = argv[2]
    fNCseed = argv[3]
    if len(argv)==5:
        jrec = int(argv[4])
    else:
        jrec = 0

    # Some strings and start/end date of Seeding input file:
    idateSeedA, idateSeedB, SeedName, SeedBatch = mjt.SeedFileTimeInfo( fNCseed, iverbose=idebug )

    # Same for model input file + time records info:
    Nt0, ztime_model, idateModA, idateModB, ModConf, ModExp = mjt.ModelFileTimeInfo( cf_uv, iverbose=idebug )

    # What records of model data can we use, based on time info from 2 input files above:
    Nt, kstrt, kstop = mjt.GetTimeSpan( rdt, ztime_model, idateSeedA, idateSeedB, idateModA, idateModB )


    cfdir = './figs/tracking'
    if iplot>0 and not path.exists(cfdir): mkdir(cfdir)
    for cd in ['nc', 'npz' ]:
        if not path.exists(cd): mkdir(cd)

    # Getting model grid metrics and friends:
    imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM = mjt.GetModelGrid( cf_mm )

    if iUVstrategy==1:
        # Get extra U,V-point metrics:
        xYv, xXv, xYu, xXu = mjt.GetModelUVGrid( cf_mm )

    # Allocating arrays for model data:
    (Nj,Ni) = np.shape( imaskt )
    xUu = np.zeros((Nj,Ni))
    xVv = np.zeros((Nj,Ni))
    xIC = np.zeros((Nj,Ni)) ; # Sea-ice concentration

    # We need a name for the intermediate backup file:
    cf_npz_itm = './npz/Initialized_buoys_'+SeedName+'.npz'


    ############################
    # Initialization / Seeding #
    ############################

    if path.exists(cf_npz_itm):
        # We save a lot of energy by using the previously generated intermediate backup file:        
        print('\n *** We found file '+cf_npz_itm+' here! So using it and skeeping first stage!')
        with np.load(cf_npz_itm) as data:
            nP    = data['nP']
            xPosG0 = data['xPosG0']
            xPosC0 = data['xPosC0']
            IDs   = data['IDs']
            vJIt  = data['vJIt']
            VRTCS = data['VRTCS']

    else:

        # Going through whole initialization / seeding process
        # ----------------------------------------------------

        with Dataset(cf_uv) as ds_UVmod:
            xIC[:,:] = ds_UVmod.variables['siconc'][0,:,:] ; # We need ice conc. at t=0 so we can cancel buoys accordingly

        zt, zIDs, XseedG, XseedC = mjt.LoadNCdataMJT( fNCseed, krec=jrec, iverbose=idebug )
        print('     => data used for seeding is read at date =',mjt.epoch2clock(zt),'\n        (shape of XseedG =',np.shape(XseedG),')')

        (nP,_) = np.shape(XseedG)

        # We want an ID for each seeded buoy:
        IDs = np.array( range(nP), dtype=int) + 1 ; # Default! No ID=0 !!!


        IDs[:] = zIDs[:]
        del zIDs

        # Find the location of each seeded buoy onto the model grid:
        nP, xPosG0, xPosC0, IDs, vJIt, VRTCS = mjt.SeedInit( IDs, XseedG, XseedC, xlatT, xlonT, xYf, xXf,
                                                             xResKM, imaskt, xIceConc=xIC, iverbose=idebug )
        del XseedG, XseedC

        # This first stage is fairly costly so saving the info:
        print('\n *** Saving intermediate data into '+cf_npz_itm+'!')
        np.savez_compressed( cf_npz_itm, nP=nP, xPosG0=xPosG0, xPosC0=xPosC0, IDs=IDs, vJIt=vJIt, VRTCS=VRTCS )

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

    if iplot>0 and idebug>1:
        mjt.ShowBuoysMap( 0, xPosG[0,:,1], xPosG[0,:,0], pvIDs=IDs, cfig=cfdir+'/INIT_Pos_buoys_'+SeedBatch+'_'+ModExp+'_'+'%4.4i'%(jt)+'.png',
                          cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1., title='IceTracker: Init Seeding' )






    ######################################
    # Loop along model data time records #
    ######################################

    # Open Input SI3 data file
    ds_UVmod = Dataset(cf_uv)

    for jt in range(Nt):

        jrec = jt + kstrt ; # access into netCDF file...

        rtmod = ds_UVmod.variables['time_counter'][jrec] ; # time of model data (center of the average period which should = rdt)
        itime = int(rtmod - rdt/2.) ; # velocitie is average under the whole rdt, at the center!
        ctime = mjt.epoch2clock(itime)
        print('\n *** Reading record #'+str(jrec)+'/'+str(Nt0)+' in SI3 file ==> date =',
              ctime,'(model:'+mjt.epoch2clock(int(rtmod))+')')
        vTime[jt] = itime

        xIC[:,:] = ds_UVmod.variables['siconc'][jrec,:,:]
        xUu[:,:] = ds_UVmod.variables['u_ice'][jrec,:,:]
        xVv[:,:] = ds_UVmod.variables['v_ice'][jrec,:,:]

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
                if   iUVstrategy == 0:
                    zU = 0.5*(xUu[jT,iT]+xUu[jT,iT-1])
                    zV = 0.5*(xVv[jT,iT]+xVv[jT-1,iT])
                    #
                elif iUVstrategy == 1:
                    # If the segment that goes from our buoy position to the F-point of the cell
                    # intersects the segment that joins the 2 V-points of the cell (v[j-1,i],v[j,i]),
                    # then it means that the nearest U-point is the one at `i-1` !
                    Fpnt = [xYf[jT,iT],xXf[jT,iT]] ; # y,x coordinates of the F-point of the cell
                    llum1 = mjt.intersect2Seg( [ry,rx], Fpnt,  [xYv[jT-1,iT],xXv[jT-1,iT]], [xYv[jT,iT],xXv[jT,iT]] )
                    llvm1 = mjt.intersect2Seg( [ry,rx], Fpnt,  [xYu[jT,iT-1],xXu[jT,iT-1]], [xYu[jT,iT],xXu[jT,iT]] )
                    if llum1:
                        zU = xUu[jT,iT-1]
                    else:
                        zU = xUu[jT,iT]
                    if llvm1:
                        zV = xVv[jT-1,iT]
                    else:
                        zV = xVv[jT,iT]
                    if idebug>1:
                        print( ' ++ Buoy position is:',ry,rx)
                        print( ' ++ position of lhs & rhs U-point:',xYu[jT,iT-1],xXu[jT,iT-1], xYu[jT,iT],xXu[jT,iT], ' llum1=',llum1)
                        print( ' ++ position of lower & upper V-point:',xYv[jT,iT-1],xXv[jT,iT-1], xYv[jT,iT],xXv[jT,iT], ' llvm1=',llvm1)

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

        print('\n')
    ### for jt in range(Nt)

    ds_UVmod.close()

    vTime[Nt] = vTime[Nt-1] + int(rdt)


    # Masking arrays:
    xPosG = np.ma.masked_where( xmask==0, xPosG )
    xPosC = np.ma.masked_where( xmask==0, xPosC )

    # ==> time to save itime, xPosXX, xPosYY, xPosLo, xPosLa into a netCDF file !
    cdt1, cdt2 = split(':',mjt.epoch2clock(vTime[0]))[0] , split(':',mjt.epoch2clock(vTime[Nt]))[0] ; # keeps at the hour precision...
    cdt1, cdt2 = str.replace( cdt1, '-', '') , str.replace( cdt2, '-', '')
    cdt1, cdt2 = str.replace( cdt1, '_', 'h') , str.replace( cdt2, '_', 'h')
    corgn = 'NEMO-SI3_'+ModConf+'_'+ModExp
    cf_nc_out = './nc/'+corgn+'_tracking_'+SeedBatch+'_'+cdt1+'_'+cdt2+'.nc'

    kk = mjt.ncSaveCloudBuoys( cf_nc_out, vTime, IDs, xPosC[:,:,0], xPosC[:,:,1], xPosG[:,:,0], xPosG[:,:,1],
                               mask=xmask[:,:,0], tunits=ctunits_expected, fillVal=FillValue, corigin=corgn )



    if iplot>0:
        # Show on the map of the Arctic:
        for jt in range(Nt+1):
            if jt%isubsamp_fig == 0:
                zLon = np.ma.masked_where( xmask[jt,:,1]==0, xPosG[jt,:,1] )
                zLat = np.ma.masked_where( xmask[jt,:,0]==0, xPosG[jt,:,0] )
                mjt.ShowBuoysMap( vTime[jt], zLon, zLat, pvIDs=IDs,
                                  cfig=cfdir+'/Pos_buoys_'+SeedBatch+'_'+ModExp+'_'+'%4.4i'%(jt)+'_'+mjt.epoch2clock(vTime[jt])+'.png',
                                  cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1.,
                                  title='IceTracker + SI3 '+ModExp+' u,v fields' )
                del zLon, zLat






    print('        => first and final dates in simulated trajectories:',mjt.epoch2clock(vTime[0]),mjt.epoch2clock(vTime[-1]),'\n')

