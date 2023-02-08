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



idebug=2

rdt = 3600. ; # time step






def debugSeeding():
    xz = np.array([
        [ 20.,84.],
        [ 50.,89.]  ])    
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


    xPosX = np.zeros((Nt,nP)) ; # x-position of buoy along the Nt records [km]
    xPosY = np.zeros((Nt,nP)) ; # y-position of buoy along the Nt records [km]
    xPosLon = np.zeros((Nt,nP))
    xPosLat = np.zeros((Nt,nP))

    lStillIn = np.zeros(nP, dtype=bool) ; # tells if a buoy is still within expected mesh/cell..

    
    vinF = np.zeros(nP, dtype=int)
    vjnF = np.zeros(nP, dtype=int)  ; # j,i indices of the F-point that defines the current mesh/cell
    #                                 # (was the nearest point when we searched for nearest point,
    #                                 #  as buoys move moves within the cell, it might not be the nearest point)

    vQuads = np.zeros(nP, dtype=Polygon) ; # stores for each buoy the polygon object associated to the current mesh/cell
    vCorners = np.zeros((nP,8), dtype=int)
    
    
    print('\n\n *** '+str(Nt)+' records in input SI3 file!\n')
    
    

    for jt in range(Nt-1):

        rtime = id_uv.variables['time_counter'][jt]
        print('\n\n *** Reading record #'+str(jt)+' in SI3 file ==> date = ', epoch2clock(int(rtime)))
        
        xUu[:,:]   = id_uv.variables['u_ice'][jt,:,:]
        xVv[:,:]   = id_uv.variables['v_ice'][jt,:,:]
    
        print('   *   current number of buoys to follow: '+str(nP))
        
        JInrstF = np.zeros((nP,2), dtype=int)
        xrjiF   = np.zeros((nP,2))
    
        for jP in range(nP):

            if jt==0:
                lStillIn[jP] = False ; #fixme: not needed!
                rlon, rlat = Xseed0G[jP,0], Xseed0G[jP,1] ; # degrees!
                rx  , ry   = Xseed0C[jP,0], Xseed0C[jP,1] ; # km !
                xPosLon[jt,jP], xPosLat[jt,jP] = rlon, rlat
                xPosX[jt,jP] , xPosY[jt,jP] = rx  , ry
            else:
                rx  , ry   = xPosX[jt,jP] , xPosY[jt,jP] ; # km !
                rlon, rlat = xPosLon[jt,jP], xPosLat[jt,jP] ; # km !

            if idebug>0:
                print('\n    * buoy (#'+str(IDs[jP])+'), jt='+str(jt)+': ry, rx =', ry, rx,'km'+': rlat, rlon =', rlat, rlon )

            if not lStillIn[jP]:
            #if 1==1:
                print('    ++++ RELOCATION NEEDED for buoy #'+str(IDs[jP])+' ++++')
                
                # We have to (re-) identify the model mesh/cell that includes current position `rlon, rlat` !
                
                # 1/ Nearest F-point on NEMO grid:
                #jnF,inF = 0,0
                #[jnF, inF] = gz.NearestPoint( (rlat,rlon), xlatF, xlonF, rd_found_km=50., j_prv=vjnF[jP], i_prv=vinF[jP] )
                [jnF, inF] = gz.NearestPoint( (rlat,rlon), xlatF, xlonF, rd_found_km=5., j_prv=0, i_prv=0 )
                vjnF[jP], vinF[jP] = jnF, inF
                #[jnF, inF] = gz.NearestPoint( (ry,rx), xYf, xXf, rd_found_km=50., j_prv=jnF, i_prv=inF )
                JInrstF[jP,:] = [ jnF, inF ]
                if idebug>1:
                    print('     ==> nearest F-point for ',rlat,rlon,' on NEMO grid:', jnF, inF, '==> lat,lon:',
                          round(xlatF[jnF,inF],3), round(xlonF[jnF,inF],3))
            
                #      o--o           x--o            o--x            o--o
                # 1 => |  | NE   2 => |  | SE    3 => |  | SW    4 => |  | NW
                #      x--o           o--o            o--o            o--x
                iq = gz.Iquadran( (rlat,rlon), xlatF, xlonF, jnF, inF, k_ew_per=-1, lforceHD=True )
                #iq = gz.Iquadran( (ry,rx), xYf, xXf, jnF, inF, k_ew_per=-1, lforceHD=True )
                if not iq in [1,2,3,4]:
                    print('PROBLEM: Fuck up for this point `iq`! `iq` = ',iq); exit(0)
                
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
                JIsSurroundMesh = gz.IDSourceMesh( (rlat,rlon), xlatF, xlonF, jnF, inF, iquadran=iq, k_ew_per=-1, lforceHD=True )
                #JIsSurroundMesh = gz.IDSourceMesh( (ry,rx), xYf, xXf, jnF, inF, iquadran=iq, k_ew_per=-1, lforceHD=True )
    
                [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = JIsSurroundMesh[:,:]
                # => indexing is anti-clockwize with (j1,i1) beeing the F-point nearest to the buoy...
    
                # The current mesh/cell as a shapely polygon object:
                MeshCell = Polygon( [ (xYf[j1,i1],xXf[j1,i1]) , (xYf[j2,i2],xXf[j2,i2]) ,
                                      (xYf[j3,i3],xXf[j3,i3]) , (xYf[j4,i4],xXf[j4,i4]) ] )                
                vQuads[jP] = MeshCell ; # store it for each buoy...
                
                if idebug>0:
                    print('     ==> vIDsrcMsh =', [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ])
                    # Buoy location:
                    if not IsInsideCell(rx, ry, MeshCell):
                        print('\nPROBLEM: buoy location is not inside the expected mesh!!!')
                        print([(xYf[j1,i1],xXf[j1,i1]), (xYf[j2,i2],xXf[j2,i2]), (xYf[j3,i3],xXf[j3,i3]), (xYf[j4,i4],xXf[j4,i4])])
                        exit(0)
                    #
                    # T-point @ center of mesh ?
                    if not IsInsideCell( xXt[jnT,inT], xYt[jnT,inT], MeshCell):
                        print('\nPROBLEM: T-point is not inside the expected mesh!!!')
                        print([(xYf[j1,i1],xXf[j1,i1]), (xYf[j2,i2],xXf[j2,i2]), (xYf[j3,i3],xXf[j3,i3]), (xYf[j4,i4],xXf[j4,i4])])
                        exit(0)
                        
                del MeshCell
                # The 4 weights for bi-linear interpolation within this cell:
                #[w1, w2, w3, w4] = gz.WeightBL( (rlat,rlon), xlatF, xlonF, JIsSurroundMesh )
                                
                # Identify the 4 corners (aka F-points) as bottom left, bottom right, upper right & and upper left
                if   iq==1:
                    vCorners[jP,:] = [ j1,i1, j2,i2, j3,i3, j4,i4 ]
                elif iq==2:
                    vCorners[jP,:] = [ j2,i2, j3,i3, j4,i4, j1,i1 ]
                elif iq==3:
                    vCorners[jP,:] = [ j3,i3, j4,i4, j1,i1, j2,i2 ]
                elif iq==4:
                    vCorners[jP,:] = [ j4,i4, j1,i1, j2,i2, j3,i3 ]


            [ jbl,ibl, jbr,ibr, jur,iur, jul,iul ] = vCorners[jP,:]                        
            print(' * lic =',jbl,ibl, jbr,ibr, jur,iur, jul,iul)
                
            # ASSUMING THAT THE ENTIRE CELL IS MOVING AT THE SAME VELOCITY: THAT OF U-POINT OF CELL
            zU, zV = xUu[jur,iur], xVv[jur,iur] ; # because the F-point is the upper-right corner
            print('    * ice velocity of the mesh: u,v =',zU, zV, 'm/s')

            # Displacement during the upcomming time step:
            dx = zU*rdt
            dy = zV*rdt
            print('      ==> displacement during `dt`: dx,dy =',dx,dy, 'm')

            # => position [km] of buoy at next time step will be:
            rx_nxt = rx + dx/1000.
            ry_nxt = ry + dy/1000.

            xPosX[jt+1,jP] = rx_nxt
            xPosY[jt+1,jP] = ry_nxt
            
            # Is it still inside our mesh:
            lStillIn[jP] = IsInsideCell(rx_nxt, ry_nxt, vQuads[jP])
            print('      ==> Still inside the same mesh???',lStillIn[jP])
            #if not lStillIn[jP]: exit(0)
            
            
        ### for jP in range(nP)

        # Updating in terms of lon,lat for all the buoys
        xPosLon[jt+1,:], xPosLat[jt+1,:] = mjt.ConvertCartesianNPSkm2Geo( xPosX[jt+1,:] , xPosY[jt+1,:] )

        #print('LOLO present and next lon, lat pos:')
        #print(xPosLon[jt,:], xPosLat[jt,:],'\n')
        #print(xPosLon[jt+1,:], xPosLat[jt+1,:],'\n')
            
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
        
