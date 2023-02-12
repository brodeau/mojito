#
import numpy as np
#
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


#from math import atan2,pi


rmin_conc = 0.6 ; # ice concentration below which we disregard the point...

rFoundKM = 5.


def IsInsideCell( py, px, CellPolygon ):
    '''
         * py, px     : global (not local to the cell) y,x position [m] or [km]
         * CellPolygon: shapely polygon object of a quadrangle mesh constructed with same unit as px and py
    '''
    pnt  = Point(py,px)
    return CellPolygon.contains(pnt)


def _ccw_( pcA, pcB, pcC ):
    '''
        * pcA: coordinates of point A => [y_A,x_A]
        * etc...
    '''
    return (pcC[0]-pcA[0])*(pcB[1]-pcA[1]) > (pcB[0]-pcA[0])*(pcC[1]-pcA[1])


def intersect2Seg( pcA, pcB, pcC, pcD ):
    '''
     Return true if line segments AB and CD intersect

        * pcA: coordinates of point A => [y_A,x_A]
        * etc...
    '''
    return ( _ccw_(pcA,pcC,pcD) != _ccw_(pcB,pcC,pcD) ) and ( _ccw_(pcA,pcB,pcC) != _ccw_(pcA,pcB,pcD) )



def Survive( kID, kjiT, pmskT, pIceC=[],  iverbose=0 ):
    '''
        Suite of tests to decide whether a buoy should be killed or not...
    '''
    (Nj,Ni) = np.shape(pmskT) ; # shape of the NEMO domain...
    
    [jT,iT] = kjiT
    ikill   = 0

    # Test: too close to NEMO domain boundaries:
    if ikill==0:
        if jT in [0,1,Nj-2,Nj-1] or iT in [0,1,Ni-2,Ni-1]:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (reaching NEMO domain boundaries)')
    
    # Test on land-sea mask / continent:
    if ikill==0:
        zmt = pmskT[jT,iT] + pmskT[jT,iT+1]+pmskT[jT+1,iT]+pmskT[jT,iT-1]+pmskT[jT-1,iT-1]
        if zmt < 5:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (too close to or over the land-sea mask)')

    # Test on sea-ice concentration:
    if ikill==0:
        if len(np.shape(pIceC))==2:
            zic = 0.2*(pIceC[jT,iT] + pIceC[jT,iT+1]+pIceC[jT+1,iT]+pIceC[jT,iT-1]+pIceC[jT-1,iT-1])
            if iverbose>1: print('     =>> 5P sea-ice concentration =',zic)
        if zic < rmin_conc:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (mean 5P ice concentration at T-point of the cell =',zic,')')
    #
    return ikill



def SeedInit( pIDs, pSG, pSC, platT, plonT, pYf, pXf, ialive, maskT, xIceConc=[], iverbose=0 ):
    '''

    '''
    from gonzag import NearestPoint
    #
    (nP,n2) = np.shape(pSG)
    #
    if np.shape(pSC)!=(nP,n2):
        print('ERROR [SeedInit]: shape disagreement for `pSG` and `pSC`!'); exit(0)
    if n2!=2:
        print('ERROR [SeedInit]: wrong shape for `pSG` and `pSC`!'); exit(0)
    #
    #                  jF,iF:       # j,i indices of the F-point that defines the current mesh/cell
    #                                 # (was the nearest point when we searched for nearest point,
    #                                 #  as buoys move moves within the cell, it might not be the nearest point)
    zjiT     = np.zeros((nP,2),   dtype=int)
    zJIvrtcs = np.zeros((nP,2,4), dtype=int)

    if np.sum(ialive) != nP:
        print('ERROR [SeedInit]: we should begin with a `ialive` filled with 1s!'); exit(0)
    
    for jP in range(nP):

        if iverbose>0:
            print('   * [SeedInit()]: focus on buoy with ID:'+str(pIDs[jP]))

        # Initial position of the buoy:
        zlat,zlon  = pSG[jP,0], pSG[jP,1] ; # degrees!
        zy  ,zx    = pSC[jP,0], pSC[jP,1] ; # km !
        zic = 1.
        # 1/ Nearest T-point on NEMO grid:
        [jT, iT] = NearestPoint( (zlat,zlon), platT, plonT, rd_found_km=rFoundKM, j_prv=0, i_prv=0 )

        if jT<0 or iT<0:
            ialive[jP] = 0
            if iverbose>0: print('        ===> I CANCEL buoy '+str(pIDs[jP])+'!!! (NO nearest T-point found for ',zlat,zlon,')')
                            
        if ialive[jP] == 1:            
            # Ok a nearest point was found!    
            if iverbose>0:
                print('     ==> nearest T-point for ',zlat,zlon,' on NEMO grid:', jT, iT, '==> lat,lon:',
                      round(platT[jT,iT],3), round(plonT[jT,iT],3))
            # Tests for canceling buoy or not:
            icncl = Survive( pIDs[jP], [jT,iT] , maskT, pIceC=xIceConc, iverbose=iverbose )
            if icncl>0: ialive[jP] = 0
            
        if ialive[jP] == 1:
            # Everything is okay, final work...
            zjiT[jP,:] = [jT,iT]
                
            # Normally, based on this nearest T-point, the point should belong
            # to the polygon defined by the 4 following F-points:
            [jf1,jf2,jf3,jf4] = [ jT-1, jT-1, jT, jT  ]
            [if1,if2,if3,if4] = [ iT-1, iT,   iT, iT-1 ]
            
            PolF = Polygon( [ (pYf[jf1,if1],pXf[jf1,if1]) , (pYf[jf2,if2],pXf[jf2,if2]) ,
                              (pYf[jf3,if3],pXf[jf3,if3]) , (pYf[jf4,if4],pXf[jf4,if4]) ] )
            lPin = IsInsideCell(zy, zx, PolF)
            if not lPin:
                print('WARNING [SeedInit()]: F-point cell not the one expected!!!'); #exit(0)
                print('        => when lookin for point:',zlat,zlon)
                ialive[jP] = 0
            else:
                # Allright => indexing is anti-clockwize, starting from bottom left F-point
                zJIvrtcs[jP,:,:] = [ [ jf1, jf2, jf3, jf4 ], [ if1, if2, if3, if4 ] ]
                        
    ### for jP in range(nP)

    return  zjiT, zJIvrtcs











def CrossedEdge( pP1, pP2, ji4vert, pY, pX,  iverbose=0 ):
    '''
        * pP1    : point 1 as [y1,x1]
        * pP2    : point 2 as [y2,x2]
        * ji4vert : array(2,4), containg the 4 j,i indices of the 4 vertices of the mesh
    '''
    #
    for kk in range(4):
        kp1 = (kk+1)%4
        [j1,i1] = ji4vert[:,kk]
        [j2,i2] = ji4vert[:,kp1]
        #
        ll = intersect2Seg( pP1, pP2,  [pY[j1,i1],pX[j1,i1]], [pY[j2,i2],pX[j2,i2]] )
        if ll: break
        #
    if iverbose>0:
        vdir = ['bottom', 'right-hand', 'upper', 'left-hand']
        print('    [CrossedEdge()]: particle is crossing the '+vdir[kk]+' edge of the mesh!')
    return kk+1




def NewHostCell( kcross, pP1, pP2, ji4vert, pY, pX,  iverbose=0 ):
    '''
       Find in which adjacent cell, the point has moved into !

       * kcross : integer describing which (old) cell edge the point has crossed

       In rare cases, the buoy could possibly move into diagonally adjacent cells
       (not only bottom,right,upper,left cells...)
       These "diagonally-adjacent" meshes yields output: 5,6,7,8
    '''
    #
    [ [ jbl, jbr, jur, jul ], [ ibl, ibr, iur, iul ] ] = ji4vert[:,:]
    #
    knhc = kcross
    if  kcross==1:
        # Crosses the bottom edge:
        if   intersect2Seg( pP1, pP2, [pY[jbl,ibl],pX[jbl,ibl]], [pY[jbl-1,ibl],pX[jbl-1,ibl]] ):
            knhc=5 ; # bottom left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jbr,ibr],pX[jbr,ibr]], [pY[jbr-1,ibr],pX[jbr-1,ibr]] ):
            knhc=6 ; # bottom right diagonal
        #
    elif kcross==2:
        # Crosses the RHS edge:
        if   intersect2Seg( pP1, pP2, [pY[jbr,ibr],pX[jbr,ibr]], [pY[jbr,ibr+1],pX[jbr,ibr+1]] ):
            knhc=6 ; # bottom right diagonal
        elif intersect2Seg( pP1, pP2, [pY[jur,iur],pX[jur,iur]], [pY[jur,iur+1],pX[jur,iur+1]] ):
            knhc=7 ; # bottom right diagonal
        #
    elif kcross==3:
        # Crosses the upper edge:
        if   intersect2Seg( pP1, pP2, [pY[jul,iul],pX[jul,iul]], [pY[jul+1,iul],pX[jul+1,iul]] ):
            knhc=8 ; # upper left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jur,iur],pX[jur,iur]], [pY[jur+1,iur],pX[jur+1,iur]] ):
            knhc=7 ; # bottom right diagonal
        #
    elif kcross==4:
        # Crosses the LHS edge:
        if   intersect2Seg( pP1, pP2, [pY[jul,iul],pX[jul,iul]], [pY[jul,iul-1],pX[jul,iul-1]] ):
            knhc=8 ; # upper left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jbl,ibl],pX[jbl,ibl]], [pY[jbl,ibl-1],pX[jbl,ibl-1]] ):
            knhc=5 ; # bottom left diagonal
    #
    if iverbose>0:
        vdir = ['bottom', 'RHS', 'upper', 'LHS', 'bottom-LHS', 'bottom-RHS', 'upper-RHS', 'upper-LHS' ]
        print('    *** Particle is moving into the '+vdir[knhc-1]+' mesh !')
    #
    return knhc



def UpdtInd4NewCell( knhc, ji4vert, kjiT ):
    '''
        Update the mesh indices according to the new host cell
    '''
    if   knhc==1:
        # went down:
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[0] = kjiT[0]-1
    elif knhc==5:
        print('LOLO: WE HAVE A 5 !!!!')
        # went left + down
        ji4vert[1,:] = ji4vert[1,:] - 1
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[1] = kjiT[1]-1
        kjiT[0] = kjiT[0]-1
    elif knhc==6:
        print('LOLO: WE HAVE A 6 !!!!')
        # went right + down
        ji4vert[1,:] = ji4vert[1,:] + 1
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[1] = kjiT[1]+1
        kjiT[0] = kjiT[0]-1
    elif knhc==2:
        # went to the right
        ji4vert[1,:] = ji4vert[1,:] + 1
        kjiT[1] = kjiT[1]+1
    elif knhc==3:
        # went up:
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[0] = kjiT[0]+1
    elif knhc==7:
        print('LOLO: WE HAVE A 7 !!!!')
        # went right + up
        ji4vert[1,:] = ji4vert[1,:] + 1
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[1] = kjiT[1]+1
        kjiT[0] = kjiT[0]+1
    elif knhc==8:
        print('LOLO: WE HAVE A 8 !!!!')
        # went right + up
        ji4vert[1,:] = ji4vert[1,:] - 1
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[1] = kjiT[1]-1
        kjiT[0] = kjiT[0]+1
    elif knhc==4:
        # went to the left
        ji4vert[1,:] = ji4vert[1,:] - 1
        kjiT[1] = kjiT[1]-1
    else:
        print('ERROR: unknown direction, knhc=',knhc)
        exit(0)

    return ji4vert, kjiT



