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






def debugSeeding():
    xz = np.array([
        [ 20.,84.]  ])
        #[ 20.,84.],
        #[ 50.,89.],
        #[100.,89.],
        #[200.,83.],
        #[300.,85.],
        # ])
    return xz
    

def IsInsideCell( px, py, CellPolygon ):
    '''
         * px, py     : global (not local to the cell) x,y position [m] or [km]
         * CellPolygon: shapely polygon object of a quadrangle mesh constructed with same unit as px and py
    '''    
    pnt  = Point(py,px)
    return CellPolygon.contains(pnt)


if __name__ == '__main__':

    if not len(argv) in [3]:
        print('Usage: '+argv[0]+' <ice_file_velocities_SI3> <mesh_mask>')
        exit(0)

    cf_uv = argv[1]
    cf_mm = argv[2]


    
    ################################################################################################
    ###                                     S E E D I N G                                        ###
    if idebug>1: Xseed0G = debugSeeding()
            
    print('\n shape of Xseed0G =',np.shape(Xseed0G))

    (nP,_) = np.shape(Xseed0G)
    IDs    = np.array( range(nP), dtype=int) + 1 ; # No ID=0 !!!

    if idebug>2:
        # Show buoys on the map at initial position (as seeded):
        mjt.ShowBuoysMap( 0, Xseed0G[:,0], Xseed0G[:,1], pvIDs=IDs, cfig='fig_initPos_seeded.png',
                          cnmfig=None, ms=10, ralpha=1., lShowDate=False )
        
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

        #xlatU  = id_mm.variables['gphiu'][0,:,:]
        #xlatV  = id_mm.variables['gphiv'][0,:,:]

        #xe1v   = id_mm.variables['e1v'][0,:,:] / 1000. ; # km !
        #xe1t   = id_mm.variables['e1t'][0,:,:] / 1000. ; # km !

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
    
    if idebug>1:
        ii = dump_2d_field( 'xXt.nc', xXt, xlon=xlonT, xlat=xlatT, name='xXt', unit='km' )
        ii = dump_2d_field( 'xXf.nc', xXf, xlon=xlonF, xlat=xlatF, name='xXf', unit='km' )

    # Xseed0C: Xseed0G converted to position in km:
    zx,zy = mjt.ConvertGeo2CartesianNPSkm(Xseed0G[:,0], Xseed0G[:,1])         
    Xseed0C = np.array([zx,zy]).T
    del zx,zy
    

    #if idebug>1:
    #    print('\n Xseed0G =',Xseed0G)
    #    zlon,zlat = mjt.ConvertCartesianNPSkm2Geo(Xseed0C[:,0], Xseed0C[:,1])
    #    ZZ = np.array([zlon,zlat]).T
    #    print('\n ZZ =',ZZ)        
    #exit(0)
    
    print('\n shape of XseedC =',np.shape(Xseed0C))
    #exit(0)
    # Recherche du pont le plus proche du pole Nord:
    #iif = np.where(xlatF>89.9)
    #print(' F-point: ',iif,'=>',xlatF[iif])
    #iit = np.where(xlatT>89.9)
    #print(' T-point: ',iit,'=>',xlatT[iit])
    #iiu = np.where(xlatU>89.98)
    #print(' U-point: ',iiu,'=>',xlatU[iiu])
    #iiv = np.where(xlatV>89.9)
    #print(' V-point: ',iiv,'=>',xlatV[iiv])
    #([jNPu],[iNPu]) = iiu
    #print(' * North Pole U-point?: xlatU[jNPu,iNPu]=',xlatU[jNPu,iNPu])
    #xXu = np.zeros((Nj,Ni))
    #xYu = np.zeros((Nj,Ni))


    
    
    #exit(0)


    

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

    
    vinF = np.zeros(nP, dtype=int)
    vjnF = np.zeros(nP, dtype=int)  ; # j,i indices of the F-point that defines the current mesh/cell
    #                                 # (was the nearest point when we searched for nearest point,
    #                                 #  as buoys move moves within the cell, it might not be the nearest point)

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
        #iq = gz.Iquadran( (ry,rx), xYf, xXf, jnF, inF, k_ew_per=-1, lforceHD=True )
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

    ilolorm=0

    #exit(0)

    for jt in range(Nt-1):

        rtime = id_uv.variables['time_counter'][jt]
        print('\n\n *** Reading record #'+str(jt)+' in SI3 file ==> date = ', epoch2clock(int(rtime)))
        
        xUu[:,:]   = id_uv.variables['u_ice'][jt,:,:]
        xVv[:,:]   = id_uv.variables['v_ice'][jt,:,:]
    
        print('   *   current number of buoys to follow: '+str(nP))
        
        for jP in range(nP):

            rx  , ry   = xPosXX[jt,jP], xPosYY[jt,jP] ; # km !
            rlon, rlat = xPosLo[jt,jP], xPosLa[jt,jP] ; # degrees !

            if idebug>0:
                print('\n    * BUOY #'+str(IDs[jP])+' => jt='+str(jt)+': ry, rx =', ry, rx,'km'+': rlat, rlon =', rlat, rlon )

            if not lStillIn[jP]:
                print('      +++ RELOCATION NEEDED for buoy #'+str(IDs[jP])+' +++')
                
                [ ibl, ibr, iur, iul ] = VERTICES_i[jP,:]
                [ jbl, jbr, jur, jul ] = VERTICES_j[jP,:]
                
                # The current mesh/cell as a shapely polygon object:
                MeshCell = Polygon( [ (xYf[jbl,ibl],xXf[jbl,ibl]) , (xYf[jbr,ibr],xXf[jbr,ibr]) ,
                                      (xYf[jur,iur],xXf[jur,iur]) , (xYf[jul,iul],xXf[jul,iul]) ] )                
                vCELLs[jP] = MeshCell ; # store it for each buoy...
                
                if idebug>0:
                    print('     ==> 4 corner points of our mesh (anti-clockwise, starting from BLC) =',
                          [ [jbl,ibl],[jbr,ibr],[jur,iur],[jul,iul] ])
                    # Buoy location:
                    if not IsInsideCell(rx, ry, MeshCell):
                        print('\nPROBLEM: buoy location is not inside the expected mesh!!!')
                        print([(xYf[jbl,ibl],xXf[jbl,ibl]), (xYf[jbr,ibr],xXf[jbr,ibr]), (xYf[jur,iur],xXf[jur,iur]),
                               (xYf[jul,iul],xXf[jul,iul])])
                        exit(0)
                    #
                    # T-point @ center of mesh ?
                    if not IsInsideCell( xXt[jnT,inT], xYt[jnT,inT], MeshCell):
                        print('\nPROBLEM: T-point is not inside the expected mesh!!!')
                        print([(xYf[jbl,ibl],xXf[jbl,ibl]), (xYf[jbr,ibr],xXf[jbr,ibr]), (xYf[jur,iur],xXf[jur,iur]),
                               (xYf[jul,iul],xXf[jul,iul])])
                        exit(0)
                        
                    print('     +++ control of location of point inside selected cell successfuly passed! :D')   
                del MeshCell

                #exit(0)
                if ilolorm==1: exit(0)
                
                # The 4 weights for bi-linear interpolation within this cell:
                #[w1, w2, w3, w4] = gz.WeightBL( (rlat,rlon), xlatF, xlonF, JIsSurroundMesh )
                                
                    
            jur = VERTICES_j[jP,2] ; # f-point is upper-right, so at 3rd position
            iur = VERTICES_j[jP,2]
            
            # ASSUMING THAT THE ENTIRE CELL IS MOVING AT THE SAME VELOCITY: THAT OF U-POINT OF CELL
            zU, zV = xUu[jur,iur], xVv[jur,iur] ; # because the F-point is the upper-right corner
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
            lSI = IsInsideCell(rx_nxt, ry_nxt, vCELLs[jP])            
            lStillIn[jP] = lSI
            print('      ==> Still inside the same mesh???',lSI)

            if not lSI:
                print('      ... need to find which wall was crossed...')
                print('   dx, dy=',dx, dy)
                [[rlon1,rlon2],[rlat1,rlat2]] = mjt.ConvertCartesianNPSkm2Geo( np.array([rx , rx_nxt]) , np.array([ry, ry_nxt]) )
                print("Present & upcomming coordinates:",rlat1,rlon1,'&',rlat2,rlon2)
                #rhdg = gz.Heading( rlat1,rlon1, rlat2,rlon2 )

                z4lat = xlatF[VERTICES_j[jP,:],VERTICES_i[jP,:]]
                z4lon = xlonF[VERTICES_j[jP,:],VERTICES_i[jP,:]]

                print(' z4lat =',z4lat)
                print(' z4lon =',z4lon)
                
                idir,cdir = gz.IQH( rlat1,rlon1, rlat2,rlon2, z4lat, z4lon )

                print('   ==> heading '+cdir+'! => idir =',idir)

                rangle_xy = atan2((ry_nxt-ry),(rx_nxt-rx))

                print('   ==> `atan2` on distances says the angle is',rangle_xy/pi,'*Pi rad')

                if rangle_xy<0. and rangle_xy>-pi/2.:
                    print('3 possible adjacent meshes: right, right+below (diag), or below! ')
                    
                    # Angles from bottom left corner mesh vertex to new location of point:
                    jrs, irs = VERTICES_j[jP,1],VERTICES_i[jP,1]                    
                    ra_rs    = atan2((ry_nxt-xYf[jrs,irs]),(rx_nxt-xYf[jrs,irs])) ; # with right side
                    print('      => angle from bottom left corner to point is',ra_rs/pi,'*Pi rad')
                    # Angles from bottom left corner mesh vertex of mesh to the right to new location of point:
                    jrsp1, irsp1 = VERTICES_j[jP,1],VERTICES_i[jP,1]                    
                    ra_rsp1    = atan2((ry_nxt-xYf[jrsp1,irsp1+1]),(rx_nxt-xYf[jrsp1,irsp1+1])) ; # with right side
                    print('      => angle from bottom left corner of next mesh to the right to point is',ra_rsp1/pi,'*Pi rad')
                    
                    
                    
                    #exit(0)
                idir=14
                
                [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = JIsSurroundMesh[jP,:,:]
                
                if   idir==1:
                    # went down:
                    VERTICES_j[jP,:] = VERTICES_j[jP,:] - 1
                elif idir==2:
                    # went to the right
                    VERTICES_i[jP,:] = VERTICES_i[jP,:] + 1
                elif idir==21:
                    # went right + down
                    VERTICES_i[jP,:] = VERTICES_i[jP,:] + 1
                    VERTICES_j[jP,:] = VERTICES_j[jP,:] - 1
                elif idir==23:
                    # went right + up
                    VERTICES_i[jP,:] = VERTICES_i[jP,:] + 1
                    VERTICES_j[jP,:] = VERTICES_j[jP,:] + 1                    
                elif idir==3:
                    # went up:
                    VERTICES_j[jP,:] = VERTICES_j[jP,:] + 1
                elif idir==4:
                    # went to the left
                    VERTICES_i[jP,:] = VERTICES_i[jP,:] - 1
                else:
                    print('ERROR: unknown direction, idir=',idir)
                    exit(0)

                ilolorm=1
                #lSIx = IsInsideCell(rx_nxt, ry,     vCELLs[jP]) ; # effect of just applying dx...
                #lSIy = IsInsideCell(rx_nxt    , ry_nxt, vCELLs[jP]) ; # effect of just applying dx...
                #print('lSIx, lSIy=',lSIx, lSIy)

                #exit(0)
            #if not lStillIn[jP]: exit(0)
            

            print('LOLO: distance to nearest wall:',vCELLs[jP].exterior.distance(Point(ry_nxt,rx_nxt)))  # 
            #exit(0)
        ### for jP in range(nP)

        # Updating in terms of lon,lat for all the buoys
        xPosLo[jt+1,:], xPosLa[jt+1,:] = mjt.ConvertCartesianNPSkm2Geo( xPosXX[jt+1,:] , xPosYY[jt+1,:] )

        #print('LOLO present and next lon, lat pos:')
        #print(xPosLo[jt,:], xPosLa[jt,:],'\n')
        #print(xPosLo[jt+1,:], xPosLa[jt+1,:],'\n')
            
        print('\n\n')



            
exit(0)
    
if idebug>3:
    #       dxU
    #     F------F
    #     !      !
    # dyL !  T.  ! dyR
    #     !      !
    #     F------F
    #       dxD
    dlonU = xlonF[jur,iur]-xlonF[jul,iul]
    dlonD = xlonF[jbr,ibr]-xlonF[jbl,ibl]
    dlatL = xlatF[jul,iul]-xlatF[jbl,ibl]
    dlatR = xlatF[jur,iur]-xlatF[jbr,ibr]    
    print(' dlonD, dlonU =',dlonD, dlonU, )
    print(' dlatL, dlatR =',dlatL, dlatR)


    print('    * buoy position known from proj on rlat,rlon: ry,rx =',ry,rx)
    ry_intrpl = w1*xYf[j1,i1] + w2*xYf[j2,i2] + w3*xYf[j3,i3] + w4*xYf[j4,i4]
    rx_intrpl = w1*xXf[j1,i1] + w2*xXf[j2,i2] + w3*xXf[j3,i3] + w4*xXf[j4,i4]
    print('    * buoy position known from interpolation from the 4 sur. F-points: ry,rx =',ry_intrpl,rx_intrpl)

#







#exit(0)


# Interpolation rather than trigonometry:
#rlat_e = w1*xlatF[j1,i1] + w2*xlatF[j2,i2] + w3*xlatF[j3,i3] + w4*xlatF[j4,i4]
#rlon_e = w1*xlonF[j1,i1] + w2*xlonF[j2,i2] + w3*xlonF[j3,i3] + w4*xlonF[j4,i4]

# In terms of jj,ji:
#xrjiF[jP,0] =  w1*float(j1) + w2*float(j2) + w3*float(j3) + w4*float(j4) ; # rjj
#xrjiF[jP,1] =  w1*float(i1) + w2*float(i2) + w3*float(i3) + w4*float(i4) ; # rji    
#print('     ==> location in terms of j,i =',xrjiF[jP,:])
        
