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

idebug=3

rdt = 3600. ; # time step

toDegrees = 180./pi

isubsamp_fig = 24 ; # frequency, in number of model records, we spawn a figure on the map (if idebug>2!!!)



#         [ 20.,84.]  ])
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
        imaskt = id_mm.variables['tmask'][0,0,:,:]
        imasku = id_mm.variables['umask'][0,0,:,:]
        imaskv = id_mm.variables['vmask'][0,0,:,:]
        xlonF  = id_mm.variables['glamf'][0,:,:]
        xlatF  = id_mm.variables['gphif'][0,:,:]
        xlonT  = id_mm.variables['glamt'][0,:,:]
        xlatT  = id_mm.variables['gphit'][0,:,:]

    (Nj,Ni) = np.shape(imaskt)

    imaskt = np.array(imaskt, dtype=int)
    imasku = np.array(imasku, dtype=int)
    imaskv = np.array(imaskv, dtype=int)

    xlonF = np.mod( xlonF, 360. )
    xlonT = np.mod( xlonT, 360. )

    xXt = np.zeros((Nj,Ni))
    xYt = np.zeros((Nj,Ni))
    xXf = np.zeros((Nj,Ni))
    xYf = np.zeros((Nj,Ni))
    xUu = np.zeros((Nj,Ni))
    xVv = np.zeros((Nj,Ni))


    # NPS projection (Cartesian, in km):
    xXt[:,:], xYt[:,:] = mjt.ConvertGeo2CartesianNPSkm(xlonT, xlatT)
    xXf[:,:], xYf[:,:] = mjt.ConvertGeo2CartesianNPSkm(xlonF, xlatF)

    if idebug>2:
        ii = dump_2d_field( 'xXt.nc', xXt, xlon=xlonT, xlat=xlatT, name='xXt', unit='km' )
        ii = dump_2d_field( 'xXf.nc', xXf, xlon=xlonF, xlat=xlatF, name='xXf', unit='km' )

    # Xseed0C: Xseed0G converted to position in km:
    zx,zy = mjt.ConvertGeo2CartesianNPSkm(Xseed0G[:,0], Xseed0G[:,1])
    Xseed0C = np.array([zx,zy]).T
    del zx,zy

    print('\n shape of XseedC =',np.shape(Xseed0C))



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

    vinF = np.zeros(nP, dtype=int)
    vjnF = np.zeros(nP, dtype=int)  ; # j,i indices of the F-point that defines the current mesh/cell
    #                                 # (was the nearest point when we searched for nearest point,
    #                                 #  as buoys move moves within the cell, it might not be the nearest point)

    vinT = np.zeros(nP, dtype=int)
    vjnT = np.zeros(nP, dtype=int)

    
    vCELLs = np.zeros(nP, dtype=Polygon) ; # stores for each buoy the polygon object associated to the current mesh/cell
    VERTICES_i = np.zeros((nP,4), dtype=int)
    VERTICES_j = np.zeros((nP,4), dtype=int)

    JIsSurroundMesh = np.zeros((nP,4,2), dtype=int)

    print('\n\n *** '+str(Nt)+' records in input SI3 file!\n')





    # Initialization, seeding:
    for jP in range(nP):

        rlon, rlat = Xseed0G[jP,0], Xseed0G[jP,1] ; # degrees!
        rx  , ry   = Xseed0C[jP,0], Xseed0C[jP,1] ; # km !
        xPosLo[0,jP], xPosLa[0,jP] = rlon, rlat
        xPosXX[0,jP], xPosYY[0,jP] =  rx ,  ry

        # 1/ Nearest F-point on NEMO grid:
        [jnF, inF] = gz.NearestPoint( (rlat,rlon), xlatF, xlonF, rd_found_km=5., j_prv=0, i_prv=0 )
        vjnF[jP], vinF[jP] = jnF, inF

        if idebug>1:
            print('     ==> nearest F-point for ',rlat,rlon,' on NEMO grid:', jnF, inF, '==> lat,lon:',
                  round(xlatF[jnF,inF],3), round(xlonF[jnF,inF],3))

        #      o--o           x--o            o--x            o--o
        # 1 => |  | NE   2 => |  | SE    3 => |  | SW    4 => |  | NW
        #      x--o           o--o            o--o            o--x
        iq = gz.Iquadran( (rlat,rlon), xlatF, xlonF, jnF, inF, k_ew_per=-1, lforceHD=True )
        if not iq in [1,2,3,4]:
            print('PROBLEM (init): Fuck up for this point `iq`! `iq` = ',iq); exit(0)

        # Indices of the T-point in the center of the mesh ():
        if   iq==1:
            jnT,inT = jnF+1,inF+1
        elif iq==2:
            jnT,inT = jnF  ,inF+1
        elif iq==3:
            jnT,inT = jnF  ,inF
        elif iq==4:
            jnT,inT = jnF+1,inF
        print('     =>> coord of center of mesh (T-point) => lat,lon =',
              round(xlatT[jnT,inT],3), round(xlonT[jnT,inT],3), '(iq =',iq,')')

        vinT[jP] = inT
        vjnT[jP] = jnT
                
        # Find the `j,i` indices of the 4 points composing the source mesh that includes the target point
        #  starting with the nearest point
        zJI4vertices = gz.IDSourceMesh( (rlat,rlon), xlatF, xlonF, jnF, inF, iquadran=iq,
                                                   k_ew_per=-1, lforceHD=True )

        [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = zJI4vertices
        # => indexing is anti-clockwize with (j1,i1) beeing the F-point nearest to the buoy...

        # Identify the 4 corners (aka F-points) as bottom left, bottom right, upper right & and upper left
        if   iq==1:
            VERTICES_i[jP,:] = [ i1, i2, i3, i4 ]
            VERTICES_j[jP,:] = [ j1, j2, j3, j4 ]
        elif iq==2:
            VERTICES_i[jP,:] = [ i2, i3, i4, i1 ]
            VERTICES_j[jP,:] = [ j2, j3, j4, j1 ]
        elif iq==3:
            VERTICES_i[jP,:] = [ i3, i4, i1, i2 ]
            VERTICES_j[jP,:] = [ j3, j4, j1, j2 ]
        elif iq==4:
            VERTICES_i[jP,:] = [ i4, i1, i2, i3 ]
            VERTICES_j[jP,:] = [ j4, j1, j2, j3 ]

    ### for jP in range(nP)



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
    
                    [ ibl, ibr, iur, iul ] = VERTICES_i[jP,:]
                    [ jbl, jbr, jur, jul ] = VERTICES_j[jP,:]
                    if idebug>0:
                        print('     ==> 4 corner points of our mesh (anti-clockwise, starting from BLC) =',
                              [ [jbl,ibl],[jbr,ibr],[jur,iur],[jul,iul] ])
    
                    # The current mesh/cell as a shapely polygon object:
                    vCELLs[jP] = Polygon( [ (xYf[jbl,ibl],xXf[jbl,ibl]) , (xYf[jbr,ibr],xXf[jbr,ibr]) ,
                                            (xYf[jur,iur],xXf[jur,iur]) , (xYf[jul,iul],xXf[jul,iul]) ] )
                    
                    if idebug>0:
                        #if idebug>1 and jt>0:
                        #    # We can have a look:
                        #    #zisrc_msh = np.array([ [VERTICES_j[jP,i],VERTICES_i[jP,i]] for i in range(4) ])
                        #    zisrc_msh = np.array([  [jbl,ibl],[jbr,ibr],[jur,iur],[jul,iul] ])
                        #    cnames    = np.array([ 'P'+str(i+1)+': '+str(VERTICES_j[jP,i])+','+str(VERTICES_i[jP,i]) for i in range(4) ], dtype='U32')
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
                    zisrc_msh = np.array([ [VERTICES_j[jP,i],VERTICES_i[jP,i]] for i in range(4) ])
                    cnames    = np.array([ 'P'+str(i+1)+': '+str(VERTICES_j[jP,i])+','+str(VERTICES_i[jP,i]) for i in range(4) ], dtype='U32')
                    #mjt.PlotMesh( (rlat,rlon), xlatF, xlonF, zisrc_msh, vnames=cnames,
                    #              fig_name='mesh_lon-lat_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                    #              pcoor_extra=(xlatT[jnT,inT],xlonT[jnT,inT]), label_extra='T-point' )
                    mjt.PlotMesh( ( ry , rx ),  xYf ,  xXf,  zisrc_msh, vnames=cnames,
                                  fig_name='mesh_X-Y_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                                  pcoor_extra=(xYt[jnT,inT],xXt[jnT,inT]), label_extra='T-point' )
    
                
                # j,i indices of the cell we are dealing with = that of the upper-right F-point!!!
                jM = VERTICES_j[jP,2]
                iM = VERTICES_i[jP,2]
    
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
                    # We mneed to find which wall was crossed...
                    icross = mjt.CrossedEdge( [ry,rx], [ry_nxt,rx_nxt], VERTICES_j[jP,:], VERTICES_i[jP,:], xYf, xXf, iverbose=idebug )
    
                    # Now, in rare cases, the buoy could possibly move into diagonally adjacent cells (not only bottom,right,upper,left cells...)
                    [ ibl, ibr, iur, iul ] = VERTICES_i[jP,:]
                    [ jbl, jbr, jur, jul ] = VERTICES_j[jP,:]
                    idir = icross
                    if  icross==1:
                        # Crosses the bottom edge:
                        if   mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jbl,ibl],xXf[jbl,ibl]], [xYf[jbl-1,ibl],xXf[jbl-1,ibl]] ):
                            idir=5 ; # bottom left diagonal
                        elif mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jbr,ibr],xXf[jbr,ibr]], [xYf[jbr-1,ibr],xXf[jbr-1,ibr]] ):
                            idir=6 ; # bottom right diagonal
                        #
                    elif icross==2:
                        # Crosses the RHS edge:
                        if   mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jbr,ibr],xXf[jbr,ibr]], [xYf[jbr,ibr+1],xXf[jbr,ibr+1]] ):
                            idir=6 ; # bottom right diagonal
                        elif mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jur,iur],xXf[jur,iur]], [xYf[jur,iur+1],xXf[jur,iur+1]] ):
                            idir=7 ; # bottom right diagonal
                        #
                    elif icross==3:
                        # Crosses the upper edge:
                        if   mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jul,iul],xXf[jul,iul]], [xYf[jul+1,iul],xXf[jul+1,iul]] ):
                            idir=8 ; # upper left diagonal
                        elif mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jur,iur],xXf[jur,iur]], [xYf[jur+1,iur],xXf[jur+1,iur]] ):
                            idir=7 ; # bottom right diagonal
                        #
                    elif icross==4:
                        # Crosses the LHS edge:
                        if   mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jul,iul],xXf[jul,iul]], [xYf[jul,iul-1],xXf[jul,iul-1]] ):
                            idir=8 ; # upper left diagonal
                        elif mjt.intersect2Seg( [ry,rx], [ry_nxt,rx_nxt], [xYf[jbl,ibl],xXf[jbl,ibl]], [xYf[jbl,ibl-1],xXf[jbl,ibl-1]] ):
                            idir=5 ; # bottom left diagonal
    
                    if idebug>0:
                        vdir = ['bottom', 'RHS', 'upper', 'LHS', 'bottom-LHS', 'bottom-RHS', 'upper-RHS', 'upper-LHS' ]
                        print('    *** Particle is moving into the '+vdir[idir-1]+' mesh !')
    
    
                    if   idir==1:
                        # went down:
                        VERTICES_j[jP,:] = VERTICES_j[jP,:] - 1
                        jnT = jnT-1
                    elif idir==5:
                        print('LOLO: WE HAVE A 5 !!!!')
                        # went left + down
                        VERTICES_i[jP,:] = VERTICES_i[jP,:] - 1
                        VERTICES_j[jP,:] = VERTICES_j[jP,:] - 1
                        inT = inT-1
                        jnT = jnT-1
                    elif idir==6:
                        print('LOLO: WE HAVE A 6 !!!!')
                        # went right + down
                        VERTICES_i[jP,:] = VERTICES_i[jP,:] + 1
                        VERTICES_j[jP,:] = VERTICES_j[jP,:] - 1
                        inT = inT+1
                        jnT = jnT-1                    
                    elif idir==2:
                        # went to the right
                        VERTICES_i[jP,:] = VERTICES_i[jP,:] + 1
                        inT = inT+1
                    elif idir==3:
                        # went up:
                        VERTICES_j[jP,:] = VERTICES_j[jP,:] + 1
                        jnT = jnT+1
                    elif idir==7:
                        print('LOLO: WE HAVE A 7 !!!!')
                        # went right + up
                        VERTICES_i[jP,:] = VERTICES_i[jP,:] + 1
                        VERTICES_j[jP,:] = VERTICES_j[jP,:] + 1
                        inT = inT+1
                        jnT = jnT+1
                    elif idir==8:
                        print('LOLO: WE HAVE A 8 !!!!')
                        # went right + up
                        VERTICES_i[jP,:] = VERTICES_i[jP,:] - 1
                        VERTICES_j[jP,:] = VERTICES_j[jP,:] + 1
                        inT = inT-1
                        jnT = jnT+1
                    elif idir==4:
                        # went to the left
                        VERTICES_i[jP,:] = VERTICES_i[jP,:] - 1
                        inT = inT-1
                    else:
                        print('ERROR: unknown direction, idir=',idir)
                        exit(0)
                    #
                    vinT[jP] = inT
                    vjnT[jP] = jnT
    
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
