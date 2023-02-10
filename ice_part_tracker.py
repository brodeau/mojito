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

idebug=2

rdt = 3600. ; # time step

toDegrees = 180./pi

isubsamp_fig = 48 ; # frequency, in number of model records, we spawn a figure on the map (if idebug>2!!!)


def debugSeeding():
    xz = np.array([
        [ 20.,84.],
        [ 50.,89.],
        [100.,85.],        
        [100.,89.],
        [180.,79.],
        [190.,75.],
        [210.,75.],        
        [200.,83.],
        [300.,85.],
                 ])
    return xz

def debugSeeding1():
    xz = np.array([
        [190.,75.],
                 ])
    return xz





if __name__ == '__main__':

    if not len(argv) in [3]:
        print('Usage: '+argv[0]+' <ice_file_velocities_SI3> <mesh_mask>')
        exit(0)

    cf_uv = argv[1]
    cf_mm = argv[2]


    ################################################################################################
    ###                                     S E E D I N G                                        ###
    if idebug in [0,1,2]: Xseed0G = debugSeeding()
    if idebug in [3]:     Xseed0G = debugSeeding1()

    print('\n shape of Xseed0G =',np.shape(Xseed0G))

    (nP,_) = np.shape(Xseed0G)
    IDs    = np.array( range(nP), dtype=int) + 1 ; # No ID=0 !!!

    ################################################################################################


    # Reading mesh metrics into mesh-mask file:
    with Dataset(cf_mm) as id_mm:
        #imaskt = id_mm.variables['tmask'][0,0,:,:]
        #imasku = id_mm.variables['umask'][0,0,:,:]
        #imaskv = id_mm.variables['vmask'][0,0,:,:]
        xlonF  = id_mm.variables['glamf'][0,:,:]
        xlatF  = id_mm.variables['gphif'][0,:,:]
        xlonT  = id_mm.variables['glamt'][0,:,:]
        xlatT  = id_mm.variables['gphit'][0,:,:]

    (Nj,Ni) = np.shape(xlonT)

    #imaskt = np.array(imaskt, dtype=int)
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


    # Conversion from Geographic coordinates (lat,lon) to Cartesian in km,
    #  ==> same North-Polar-Stereographic projection as RGPS data...
    xXt[:,:], xYt[:,:] = mjt.ConvertGeo2CartesianNPSkm(xlonT, xlatT)
    xXf[:,:], xYf[:,:] = mjt.ConvertGeo2CartesianNPSkm(xlonF, xlatF)
    
    # same for seeded initial positions, Xseed0G->Xseed0C:
    zx,zy = mjt.ConvertGeo2CartesianNPSkm(Xseed0G[:,0], Xseed0G[:,1])
    Xseed0C = np.array([zx,zy]).T
    del zx,zy

    #if idebug>2:
    #    ii = dump_2d_field( 'xXt.nc', xXt, xlon=xlonT, xlat=xlatT, name='xXt', unit='km' )
    #    ii = dump_2d_field( 'xXf.nc', xXf, xlon=xlonF, xlat=xlatF, name='xXf', unit='km' )



    #######################
    # Input SI3 data file #
    #######################

    id_uv = Dataset(cf_uv)

    Nt = id_uv.dimensions['time_counter'].size

    xPosXX = np.zeros((Nt,nP)) ; # x-position of buoy along the Nt records [km]
    xPosYY = np.zeros((Nt,nP)) ; # y-position of buoy along the Nt records [km]
    xPosLo = np.zeros((Nt,nP))
    xPosLa = np.zeros((Nt,nP))

    lStillIn = np.zeros(nP, dtype=bool) ; # tells if a buoy is still within expected mesh/cell..
    IsAlive  = np.zeros(nP, dtype=int)+1 ; # tells if a buoy is alive (1) or zombie (0) (discontinued)
    
    vCELLs = np.zeros(nP, dtype=Polygon) ; # stores for each buoy the polygon object associated to the current mesh/cell

    JIsSurroundMesh = np.zeros((nP,4,2), dtype=int)

    print('\n\n *** '+str(Nt)+' records in input SI3 file!\n')



    # Initialization, seeding:
    vjnT, vinT, VRTCS_j, VRTCS_i = mjt.SeedInit( Xseed0G, Xseed0C, xlatF, xlonF, xlatT, xlonT, iverbose=idebug )

    xPosLo[0,:], xPosLa[0,:] = Xseed0G[:,0], Xseed0G[:,1] ; # degrees!
    xPosXX[0,:], xPosYY[0,:] = Xseed0C[:,0], Xseed0C[:,1] ; # km !


    ######################################
    # Loop along model data time records #
    ######################################

    for jt in range(Nt-1):
        rtmod = id_uv.variables['time_counter'][jt] ; # time of model data (center of the average period which should = rdt)
        itime = int(rtmod - rdt/2.) ; # velocitie is average under the whole rdt, at the center!
        print('\n *** Reading record #'+str(jt)+'/'+str(Nt-1)+' in SI3 file ==> date =',
              epoch2clock(itime),'(model:'+epoch2clock(int(rtmod))+')')

        xUu[:,:]   = id_uv.variables['u_ice'][jt,:,:]
        xVv[:,:]   = id_uv.variables['v_ice'][jt,:,:]

        print('   *   current number of buoys to follow: '+str(nP))
        
        for jP in range(nP):

            if IsAlive[jP]==1:
            
                rx  , ry   = xPosXX[jt,jP], xPosYY[jt,jP] ; # km !
                rlon, rlat = xPosLo[jt,jP], xPosLa[jt,jP] ; # degrees !
    
                inT, jnT = vinT[jP] , vjnT[jP] 
                
                if idebug>0:
                    print('\n    * BUOY #'+str(IDs[jP])+' => jt='+str(jt)+': ry, rx =', ry, rx,'km'+': rlat, rlon =', rlat, rlon )
    
                ######################### N E W   M E S H   R E L O C A T I O N #################################
                if not lStillIn[jP]:
                    print('      +++ RELOCATION NEEDED for buoy #'+str(IDs[jP])+' +++')
    
                    [ ibl, ibr, iur, iul ] = VRTCS_i[jP,:]
                    [ jbl, jbr, jur, jul ] = VRTCS_j[jP,:]
                    if idebug>0:
                        print('     ==> 4 corner points of our mesh (anti-clockwise, starting from BLC) =',
                              [ [jbl,ibl],[jbr,ibr],[jur,iur],[jul,iul] ])
    
                    # The current mesh/cell as a shapely polygon object:
                    vCELLs[jP] = Polygon( [ (xYf[jbl,ibl],xXf[jbl,ibl]) , (xYf[jbr,ibr],xXf[jbr,ibr]) ,
                                            (xYf[jur,iur],xXf[jur,iur]) , (xYf[jul,iul],xXf[jul,iul]) ] )
                    
                    if idebug>0:
                        #if idebug>1 and jt>0:
                        #    # We can have a look:
                        #    #zisrc_msh = np.array([ [VRTCS_j[jP,i],VRTCS_i[jP,i]] for i in range(4) ])
                        #    zisrc_msh = np.array([  [jbl,ibl],[jbr,ibr],[jur,iur],[jul,iul] ])
                        #    cnames    = np.array([ 'P'+str(i+1)+': '+str(VRTCS_j[jP,i])+','+str(VRTCS_i[jP,i]) for i in range(4) ], dtype='U32')
                        #    mjt.PlotMesh( (rlat,rlon), xlatF, xlonF, zisrc_msh, vnames=cnames,
                        #                  fig_name='zmesh_lon-lat_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                        #                  pcoor_extra=(xlatT[jnT,inT],xlonT[jnT,inT]), label_extra='T-point' )
                        #    mjt.PlotMesh( ( ry , rx ),  xYf ,  xXf,  zisrc_msh, vnames=cnames,
                        #                  fig_name='zmesh_X-Y_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                        #                  pcoor_extra=(xYt[jnT,inT],xXt[jnT,inT]), label_extra='T-point' )
                        #
                        # Buoy location:
                        if not mjt.IsInsideCell(rx, ry, vCELLs[jP]):
                            print('\nPROBLEM: buoy location is not inside the expected mesh!!!')
                            print([(xYf[jbl,ibl],xXf[jbl,ibl]), (xYf[jbr,ibr],xXf[jbr,ibr]), (xYf[jur,iur],xXf[jur,iur]),
                                   (xYf[jul,iul],xXf[jul,iul])])
                            #print('==> distance to nearest wall:',vCELLs[jP].exterior.distance(Point(ry,rx)))
                            exit(0)
                        #
                        # T-point @ center of mesh ?
                        if not mjt.IsInsideCell( xXt[jnT,inT], xYt[jnT,inT], vCELLs[jP]):
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
                    # We can have a look:
                    zisrc_msh = np.array([ [VRTCS_j[jP,i],VRTCS_i[jP,i]] for i in range(4) ])
                    cnames    = np.array([ 'P'+str(i+1)+': '+str(VRTCS_j[jP,i])+','+str(VRTCS_i[jP,i]) for i in range(4) ], dtype='U32')
                    #mjt.PlotMesh( (rlat,rlon), xlatF, xlonF, zisrc_msh, vnames=cnames,
                    #              fig_name='mesh_lon-lat_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                    #              pcoor_extra=(xlatT[jnT,inT],xlonT[jnT,inT]), label_extra='T-point' )
                    mjt.PlotMesh( ( ry , rx ),  xYf ,  xXf,  zisrc_msh, vnames=cnames,
                                  fig_name='mesh_X-Y_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                                  pcoor_extra=(xYt[jnT,inT],xXt[jnT,inT]), label_extra='T-point' )
    
                
                # j,i indices of the cell we are dealing with = that of the upper-right F-point!!!
                jM = VRTCS_j[jP,2]
                iM = VRTCS_i[jP,2]
    
                # ASSUMING THAT THE ENTIRE CELL IS MOVING AT THE SAME VELOCITY: THAT OF U-POINT OF CELL
                zU, zV = xUu[jM,iM], xVv[jM,iM] ; # because the F-point is the upper-right corner
                if idebug>0:
                    print('    =>> read velocity at ji,jj=',iM,jM)
                    print('    * ice velocity of the mesh: u,v =',zU, zV, 'm/s')
    
                # Displacement during the upcomming time step:
                dx = zU*rdt
                dy = zV*rdt
                print('      ==> displacement during `dt`: dx,dy =',dx,dy, 'm')
    
                # => position [km] of buoy after this time step will be:
                rx_nxt = rx + dx/1000. ; # [km]
                ry_nxt = ry + dy/1000. ; # [km]
    
                xPosXX[jt+1,jP] = rx_nxt
                xPosYY[jt+1,jP] = ry_nxt
    
                # Is it still inside our mesh:
                lSI = mjt.IsInsideCell(rx_nxt, ry_nxt, vCELLs[jP])
                lStillIn[jP] = lSI
                if idebug>0: print('      ==> Still inside the same mesh???',lSI)
    
                if not lSI:                    
                    # Tells which of the 4 cell walls the point has crossed:
                    icross = mjt.CrossedEdge( [ry,rx], [ry_nxt,rx_nxt], VRTCS_j[jP,:], VRTCS_i[jP,:], xYf, xXf, iverbose=idebug )

                    # Tells in which adjacent cell the point has moved:
                    inhc = mjt.NewHostCell( icross, [ry,rx], [ry_nxt,rx_nxt], VRTCS_j[jP,:], VRTCS_i[jP,:], xYf, xXf,  iverbose=idebug )

                    # Update the mesh indices according to the new host cell:
                    VRTCS_j[jP,:],VRTCS_i[jP,:],vjnT[jP],vinT[jP] = mjt.UpdtInd4NewCell( inhc, VRTCS_j[jP,:], VRTCS_i[jP,:], jnT, inT )

                ### if not lSI

            ### if IsAlive[jP]==1
                
        ### for jP in range(nP)

        # Updating in terms of lon,lat for all the buoys at once:
        xPosLo[jt+1,:], xPosLa[jt+1,:] = mjt.ConvertCartesianNPSkm2Geo( xPosXX[jt+1,:] , xPosYY[jt+1,:] )

        if idebug>1 and jt%isubsamp_fig==0:
            # Show buoys on the map:
            mjt.ShowBuoysMap( itime,  xPosLo[jt,:], xPosLa[jt,:], pvIDs=IDs, cfig='Pos_buoys_'+'%4.4i'%(jt)+'.png',
                              cnmfig=None, ms=15, ralpha=1., lShowDate=True, zoom=1.2 )
            
        print('\n\n')
        
    ### for jt in range(Nt-1)


    id_uv.close()




#print('LOLO: distance to nearest wall:',vCELLs[jP].exterior.distance(Point(ry_nxt,rx_nxt)))  #
