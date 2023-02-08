#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

#TODO:


from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

from netCDF4 import Dataset

import gonzag as gz

from climporn import dump_2d_field, epoch2clock, clock2epoch
import mojito   as mjt

idebug=2


if idebug>0:
    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon




def debugSeeding():
    xz = np.array([
        [ 20.,84.],
        [ 50.,89.],    
        [100.,89.],
        [200.,83.],
        [300.,85.],
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
    if idebug>1: Xseed0 = debugSeeding()
            
    print('\n shape of Xseed0 =',np.shape(Xseed0))

    (nP,_) = np.shape(Xseed0)
    IDs    = np.array( range(nP), dtype=int) + 1 ; # No ID=0 !!!

    if idebug>2:
        # Show buoys on the map at initial position (as seeded):
        mjt.ShowBuoysMap( 0, Xseed0[:,0], Xseed0[:,1], pvIDs=IDs, cfig='fig_initPos_seeded.png',
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


    # NPS projection (Cartesian, in km):
    xXt[:,:], xYt[:,:] = mjt.ConvertGeotoCartesianNPS(xlonT, xlatT)
    xXf[:,:], xYf[:,:] = mjt.ConvertGeotoCartesianNPS(xlonF, xlatF)
    
    if idebug>1:
        ii = dump_2d_field( 'xXt.nc', xXt, xlon=xlonT, xlat=xlatT, name='xXt', unit='km' )
        ii = dump_2d_field( 'xXf.nc', xXf, xlon=xlonF, xlat=xlatF, name='xXf', unit='km' )

    # Xseed0 
    zx,zy = mjt.ConvertGeotoCartesianNPS(Xseed0[:,0], Xseed0[:,1]) 
        
    Xseed0C = np.array([zx,zy]).T

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

    print('\n\n *** '+str(Nt)+' records in input SI3 file!')
    
    

    for jt in range(Nt):

        print('\n *** Reading record #'+str(jt)+' in SI3 file...')

        rtime = id_uv.variables['time_counter'][jt]

        print('     ==> date = ', epoch2clock(int(rtime)))
        
        xUu   = id_uv.variables['u_ice'][jt,:,:]
        xVv   = id_uv.variables['v_ice'][jt,:,:]
    
        print(' * We have '+str(nP)+' seeded buoys to follow!')
        
        JInrstF = np.zeros((nP,2), dtype=int)
        xrjiF   = np.zeros((nP,2))
    
        for jP in range(nP):
    
            rlon, rlat =  Xseed0[jP,0],  Xseed0[jP,1] ; # degrees!
            rx  , ry   = Xseed0C[jP,0], Xseed0C[jP,1] ; # km !

            if idebug>0: print('\n *** New buoy (#'+str(jP)+'): lat, lon =', rlat, rlon )
    
    
            # 1/ find nearest F-point on NEMO grid:
            jnF = 0
            inF = 0
            [jnF, inF] = gz.NearestPoint( (rlat,rlon), xlatF, xlonF, rd_found_km=50., j_prv=jnF, i_prv=inF )
            JInrstF[jP,:] = [ jnF, inF ]
            #
            if idebug>1:
                print('     ==> nearest point on NEMO grid:', jnF, inF)
                print('          ==> double check:', xlatF[jnF,inF], xlonF[jnF,inF])
    
    
            #      o--o           x--o            o--x            o--o
            # 1 => |  | NE   2 => |  | SE    3 => |  | SW    4 => |  | NW
            #      x--o           o--o            o--o            o--x
            iq = gz.Iquadran( (rlat,rlon), xlatF, xlonF, jnF, inF, k_ew_per=-1, lforceHD=True )
            print(' *** iq = ',iq)


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
            print('  =>> coord of center of mesh (T-point) => lat,lon =',xlatT[jnT,inT], xlonT[jnT,inT])
            
            # Find the `j,i` indices of the 4 points composing the source mesh that includes the target point
            #  starting with the nearest point
    
            JIsSurroundMesh = gz.IDSourceMesh( (rlat,rlon), xlatF, xlonF, jnF, inF, iquadran=iq, k_ew_per=-1, lforceHD=True )
            
            [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = JIsSurroundMesh[:,:]
    
            if idebug>0: print('     ==> vIDsrcMsh =\n', [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ])
    

            
            if idebug>0:
                
                MeshCell = Polygon( [ (xYf[j1,i1],xXf[j1,i1]) , (xYf[j2,i2],xXf[j2,i2]) ,
                                      (xYf[j3,i3],xXf[j3,i3]) , (xYf[j4,i4],xXf[j4,i4]) ] )

                # Buoy location:
                pointB   = Point(ry,rx)
                if not MeshCell.contains(pointB):
                    print('\nPROBLEM: point `rx,ry` is not inside the expected mesh!!!')
                    print([(xYf[j1,i1],xXf[j1,i1]), (xYf[j2,i2],xXf[j2,i2]), (xYf[j3,i3],xXf[j3,i3]), (xYf[j4,i4],xXf[j4,i4])])
                    exit(0)
                #
                # T-point @ center of mesh:
                ryT,rxT = xYt[jnT,inT],xXt[jnT,inT]
                pointT   = Point(ryT,rxT)
                if not MeshCell.contains(pointT):
                    print('\nPROBLEM: T-point is not inside the expected mesh!!!')
                    print([(xYf[j1,i1],xXf[j1,i1]), (xYf[j2,i2],xXf[j2,i2]), (xYf[j3,i3],xXf[j3,i3]), (xYf[j4,i4],xXf[j4,i4])])
                    exit(0)
                #
                del pointB, pointT



            exit(0)


            [w1, w2, w3, w4] = gz.WeightBL( (rlat,rlon), xlatF, xlonF, JIsSurroundMesh )
    






            exit(0)

            
            # Interpolation rather than trigonometry:
            #rlat_e = w1*xlatF[j1,i1] + w2*xlatF[j2,i2] + w3*xlatF[j3,i3] + w4*xlatF[j4,i4]
            #rlon_e = w1*xlonF[j1,i1] + w2*xlonF[j2,i2] + w3*xlonF[j3,i3] + w4*xlonF[j4,i4]
    
            # In terms of jj,ji:
            xrjiF[jP,0] =  w1*float(j1) + w2*float(j2) + w3*float(j3) + w4*float(j4) ; # rjj
            xrjiF[jP,1] =  w1*float(i1) + w2*float(i2) + w3*float(i3) + w4*float(i4) ; # rji
    
            print('     ==> location in terms of j,i =',xrjiF[jP,:])
        
