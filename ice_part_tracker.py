#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

#TODO:
#     * gz.NearestPoint should be used only at seeding time / jt==0
#     Once we detect that a move makes the buoy leave the current mesh
#     we should be able to give a number between 1 and 4 to tell
#     through which cell wall the crossing happen, thus it should be
#     easy in upcomming itteration to update the indices for the new cell
#     that contains the buoy !!!


from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

from netCDF4 import Dataset

import gonzag as gz

from climporn import dump_2d_field, epoch2clock, clock2epoch
import mojito   as mjt

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from math import atan2,pi

idebug=0
iplot=1

rdt = 3600. ; # time step

toDegrees = 180./pi

isubsamp_fig = 48 ; # frequency, in number of model records, we spawn a figure on the map (if idebug>2!!!)

#seeding_type='nemo_Tpoint' ; iHSS=20
#seeding_type='debug'
seeding_type='mojitoNC' ; fNCseed = '/home/laurent/DEV/mojito/TEST_rgps/RGPS_ice_drift_1997-01-01_1997-02-01_lb.nc'

ctunits_expected = 'seconds since 1970-01-01 00:00:00' ; # we expect UNIX/EPOCH time in netCDF files!


def debugSeeding():
    zLatLon = np.array([
        [84. ,  20.],
        [89. ,  50.],
        [63. , -11.], # point outside of domain
        [85. , 100.],        
        [89. , 100.],
        [79. , 180.],
        [76. ,  46.], # No ice, Barents Sea
        [75. , 190.],
        [85.2, -15.], # Too close to land
        [75. , 210.],
        [75. , -72.], # Baffin bay
        [83. , 200.],
        [79. , -42.], # over mask (Greenland!)
        [85. , 300.]  ])
    return zLatLon

def debugSeeding1():
    zLatLon = np.array([ [75.,190.] ])
    return zLatLon


def nemoSeed( pmskT, platT, plonT, pIC, khss=1, fmsk_rstrct=None ):

    zmsk = pmskT[::khss,::khss]
    zlat = platT[::khss,::khss]
    zlon = plonT[::khss,::khss]

    (Nj,Ni) = np.shape(zmsk)
    ztmp = np.zeros((Nj,Ni))
    
    #if fmsk_rstrct:
    #    with Dataset(fmsk_rstrct) as id_mr:
    #        maskR = id_mr.variables['mask'][::khss,::khss]
    #        if np.shape(maskR) != (Nj,Ni):
    #            print('ERROR [nemoSeed()]: restricted area mask does not agree in shape with model output!'); exit(0)
    
    msk_T = np.zeros((Nj,Ni), dtype='i1')

    msk_T[:,:] = zmsk[:,:]

    # Only north of ...
    msk_T[np.where(zlat < 55.)] = 0
        
    # Only over a decent concentration of ice:
    ztmp[:,:] = pIC[::khss,::khss]
    msk_T[np.where(ztmp < 0.9)] = 0

    (idy_keep, idx_keep) = np.where( msk_T==1 )
    NbIPt = len(idy_keep)
    zLatLon = np.zeros((NbIPt,2))
    
    for jp in range(NbIPt):
        jj = idy_keep[jp]
        ji = idx_keep[jp]
        zLatLon[jp,:] = [ zlat[jj,ji],zlon[jj,ji] ]
        
    del zmsk,zlat,zlon
        
    return zLatLon


if __name__ == '__main__':

    if not len(argv) in [3]:
        print('Usage: '+argv[0]+' <ice_file_velocities_SI3> <mesh_mask>')
        exit(0)

    
        
    cf_uv = argv[1]
    cf_mm = argv[2]




    # Reading mesh metrics into mesh-mask file:
    with Dataset(cf_mm) as id_mm:
        imaskt = id_mm.variables['tmask'][0,0,:,:]
        #imasku = id_mm.variables['umask'][0,0,:,:]
        #imaskv = id_mm.variables['vmask'][0,0,:,:]
        xlonF  = id_mm.variables['glamf'][0,:,:]
        xlatF  = id_mm.variables['gphif'][0,:,:]
        xlonT  = id_mm.variables['glamt'][0,:,:]
        xlatT  = id_mm.variables['gphit'][0,:,:]

    (Nj,Ni) = np.shape(xlonT)

    imaskt = np.array(imaskt, dtype=int)
    #imasku = np.array(imasku, dtype=int)
    #imaskv = np.array(imaskv, dtype=int)

    xlonF = np.mod( xlonF, 360. )
    xlonT = np.mod( xlonT, 360. )

    xXt = np.zeros((Nj,Ni))
    xYt = np.zeros((Nj,Ni))
    xXf = np.zeros((Nj,Ni))
    xYf = np.zeros((Nj,Ni))
    xUu = np.zeros((Nj,Ni))
    xVv = np.zeros((Nj,Ni))
    xIC = np.zeros((Nj,Ni)) ; # Sea-ice concentration

    # Conversion from Geographic coordinates (lat,lon) to Cartesian in km,
    #  ==> same North-Polar-Stereographic projection as RGPS data...
    xYt[:,:], xXt[:,:] = mjt.ConvertGeo2CartesianNPSkm(xlatT, xlonT)
    xYf[:,:], xXf[:,:] = mjt.ConvertGeo2CartesianNPSkm(xlatF, xlonF)
    del xlatF, xlonF




    

    #if idebug>2:
    #    ii = dump_2d_field( 'xXt.nc', xXt, xlon=xlonT, xlat=xlatT, name='xXt', unit='km' )
    #    ii = dump_2d_field( 'xXf.nc', xXf, xlon=xlonF, xlat=xlatF, name='xXf', unit='km' )



    #######################
    # Input SI3 data file #
    #######################

    id_uv = Dataset(cf_uv)

    Nt = id_uv.dimensions['time_counter'].size

    print('\n\n *** '+str(Nt)+' records in input SI3 file!\n')


    






    ############################
    # Initialization / Seeding #
    ############################
    
    # We need the sea-ice concentration at t=0 so we can cancel buoys if no ice:
    xIC[:,:] = id_uv.variables['siconc'][0,:,:]


    ################################################################################################
    ###                                     S E E D I N G                                        ###
    if seeding_type=='nemo_Tpoint':
        Xseed0G = nemoSeed( imaskt, xlatT, xlonT, xIC, khss=iHSS, fmsk_rstrct=None )
        #
    elif seeding_type=='debug':
        if idebug in [0,1,2]: Xseed0G = debugSeeding()
        if idebug in [3]:     Xseed0G = debugSeeding1()
        #
    elif seeding_type=='mojitoNC':
        zt, zIDs, zlat, zlon, zy, zx = mjt.LoadDataMJT( fNCseed, krec=0, iverbose=idebug )
        Xseed0G = np.array([zlat,zlon]).T
        Xseed0C = np.array([zy,zx]).T
        print(' *** Read seeding positions in '+fNCseed+' !')
        print('     => date =',epoch2clock(zt))
        #print(np.shape(Xseed0G), np.shape(Xseed0G))
        del zt, zlat, zlon, zy, zx
        #print('Xseed0G =',Xseed0G)
        #print('Xseed0C =',Xseed0C)
        #exit(0)
        
    print('\n shape of Xseed0G =',np.shape(Xseed0G))

    (nP,_) = np.shape(Xseed0G)

    if seeding_type!='mojitoNC':
        # same for seeded initial positions, Xseed0G->Xseed0C:
        zy,zx = mjt.ConvertGeo2CartesianNPSkm(Xseed0G[:,0], Xseed0G[:,1])
        Xseed0C = np.array([zy,zx]).T
        del zx,zy

    ################################################################################################



    
    # Vectors with initial number of buoys
    IDs      = np.array( range(nP), dtype=int) + 1 ; # No ID=0 !!!
    IsAlive  = np.zeros(       nP , dtype=int) + 1 ; # tells if a buoy is alive (1) or zombie (0) (discontinued)

    if seeding_type=='mojitoNC':
        IDs[:] = zIDs[:]
        del zIDs
    
    vJInT, VRTCS = mjt.SeedInit( IDs, Xseed0G, Xseed0C, xlatT, xlonT, xYf, xXf, IsAlive, imaskt,
                                 xIceConc=xIC, iverbose=idebug )

    # Getting rid of useless seeded buoys
    nPn = np.sum(IsAlive)
    if nPn<nP:
        (idxGone,) = np.where(IsAlive==0)
        print('\n *** '+str(nP-nPn)+' "to-be-seeded" "buoys have been canceled by `SeedInit()`')
        print('        => their IDs:',IDs[idxGone])
        print('      ==> need to shrink some arrays, number of valid buoys is now',nPn)
        (idxKeep,) = np.where(IsAlive==1)
        IDs     = IDs[idxKeep]
        IsAlive = np.zeros(nP,dtype=int) + 1
        vJInT   = vJInT[idxKeep,:]
        VRTCS   = VRTCS[idxKeep,:,:]
        nP = nPn
        del nPn,idxGone
    else:
        idxKeep = np.arange(nP,dtype=int)

    # Allocation for nP buoys:
    vTime  = np.zeros( Nt+1, dtype=int ) ; # UNIX epoch time associated to position below
    xmaskB = np.zeros((Nt+1,nP), dtype='i1')
    xPosXX = np.zeros((Nt+1,nP)) -9999.  ; # x-position of buoy along the Nt records [km]
    xPosYY = np.zeros((Nt+1,nP)) -9999.  ; # y-position of buoy along the Nt records [km]
    xPosLo = np.zeros((Nt+1,nP)) -9999.
    xPosLa = np.zeros((Nt+1,nP)) -9999.
    
    vCELLs   = np.zeros(nP, dtype=Polygon) ; # stores for each buoy the polygon object associated to the current mesh/cell
    lStillIn = np.zeros(nP, dtype=bool) ; # tells if a buoy is still within expected mesh/cell..
    #-----------------------------------------------------------
    
    
    xPosLa[0,:], xPosLo[0,:] = Xseed0G[idxKeep,0], Xseed0G[idxKeep,1] ; # degrees!
    xPosYY[0,:], xPosXX[0,:] = Xseed0C[idxKeep,0], Xseed0C[idxKeep,1] ; # km !
    xmaskB[0,:] = 1
    del Xseed0G, Xseed0C, idxKeep

    if iplot>0:
        if not path.exists('./figs'): mkdir('./figs')
        mjt.ShowBuoysMap( 0,  xPosLo[0,:], xPosLa[0,:], pvIDs=IDs, cfig='./figs/INIT_Pos_buoys_'+'%4.4i'%(0)+'.png',
                          cnmfig=None, ms=15, ralpha=1., lShowDate=False, zoom=1.2 )

    
    ######################################
    # Loop along model data time records #
    ######################################

    for jt in range(Nt):
        rtmod = id_uv.variables['time_counter'][jt] ; # time of model data (center of the average period which should = rdt)
        itime = int(rtmod - rdt/2.) ; # velocitie is average under the whole rdt, at the center!
        print('\n *** Reading record #'+str(jt)+'/'+str(Nt-1)+' in SI3 file ==> date =',
              epoch2clock(itime),'(model:'+epoch2clock(int(rtmod))+')')
        vTime[jt] = itime

        
        xIC[:,:] = id_uv.variables['siconc'][jt,:,:]
        xUu[:,:] = id_uv.variables['u_ice'][jt,:,:]
        xVv[:,:] = id_uv.variables['v_ice'][jt,:,:]

        print('   *   current number of buoys to follow: '+str(nP))
        
        for jP in range(nP):

            if IsAlive[jP]==1:
            
                rx  , ry   = xPosXX[jt,jP], xPosYY[jt,jP] ; # km !
                rlon, rlat = xPosLo[jt,jP], xPosLa[jt,jP] ; # degrees !
    
                [jnT,inT] = vJInT[jP,:]; # indices for T-point at center of current cell
                
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
    
                # The 4 weights for bi-linear interpolation within this cell:
                #[w1, w2, w3, w4] = gz.WeightBL( (rlat,rlon), xlatF, xlonF, JIsSurroundMesh )

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
                [jT,iT] = vJInT[jP,:]
                
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
    
                xPosXX[jt+1,jP] = rx_nxt
                xPosYY[jt+1,jP] = ry_nxt
                xmaskB[jt+1,jP] = 1
                
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
                    VRTCS[jP,:,:],vJInT[jP,:] = mjt.UpdtInd4NewCell( inhc, VRTCS[jP,:,:], vJInT[jP,:] )

                    # Based on updated new indices, some buoys might get killed:
                    icncl = mjt.Survive( IDs[jP], vJInT[jP,:], imaskt, pIceC=xIC,  iverbose=idebug )
                    if icncl>0: IsAlive[jP]=0
                    
                ### if not lSI

            ### if IsAlive[jP]==1
                
        ### for jP in range(nP)

        # Updating in terms of lon,lat for all the buoys at once:
        xPosLa[jt+1,:], xPosLo[jt+1,:] = mjt.ConvertCartesianNPSkm2Geo( xPosYY[jt+1,:], xPosXX[jt+1,:] )

        print('\n\n')
        
    ### for jt in range(Nt-1)

    id_uv.close()

    vTime[Nt] = vTime[Nt-1] + int(rdt)


    # Masking arrays:
    xPosLo = np.ma.masked_where( xmaskB==0, xPosLo )
    xPosLa = np.ma.masked_where( xmaskB==0, xPosLa )
    xPosXX = np.ma.masked_where( xmaskB==0, xPosXX )
    xPosYY = np.ma.masked_where( xmaskB==0, xPosYY )
    
    # ==> time to save itime, xPosXX, xPosYY, xPosLo, xPosLa into a netCDF file !
    kk = mjt.ncSaveCloudBuoys( 'ice_tracking.nc', vTime, IDs, xPosYY, xPosXX, xPosLa, xPosLo, tunits=ctunits_expected, fillVal=-9999. )

    if iplot>0:
        # Show on the map of the Arctic:
        for jt in range(0,Nt,isubsamp_fig):
            mjt.ShowBuoysMap( vTime[jt],  xPosLo[jt,:], xPosLa[jt,:], pvIDs=IDs, cfig='./figs/Pos_buoys_'+'%4.4i'%(jt)+'.png',
                              cnmfig=None, ms=15, ralpha=1., lShowDate=True, zoom=1.2 )
            
