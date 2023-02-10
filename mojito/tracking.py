#
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#from math import atan2,pi


rmin_conc = 0.6 ; # ice concentration below which we disregard the point...

rFoundKM = 5.

def SeedInit( pSG, pSC, platF, plonF, platT, plonT, ialive, maskT, xIceConc=[], iverbose=0 ):
    '''

    '''
    import numpy as np
    from gonzag import NearestPoint, Iquadran, IDSourceMesh
    #
    (nP,n2) = np.shape(pSG)
    #
    if np.shape(pSC)!=(nP,n2):
        print('ERROR [SeedInit]: shape disagreement for `pSG` and `pSC`!'); exit(0)
    if n2!=2:
        print('ERROR [SeedInit]: wrong shape for `pSG` and `pSC`!'); exit(0)
    #
    #                  jnF,inF:       # j,i indices of the F-point that defines the current mesh/cell
    #                                 # (was the nearest point when we searched for nearest point,
    #                                 #  as buoys move moves within the cell, it might not be the nearest point)
    zjinT    = np.zeros((nP,2),   dtype=int)
    zJIvrtcs = np.zeros((nP,2,4), dtype=int)

    for jP in range(nP):

        if iverbose>1:
            print('   * [SeedInit()]: focus on buoy #'+str(jP))
        
        zlon, zlat = pSG[jP,0], pSG[jP,1] ; # degrees!
        zx  , zy   = pSC[jP,0], pSC[jP,1] ; # km !
        zic = 1.
        # 1/ Nearest F-point on NEMO grid:
        [jnF, inF] = NearestPoint( (zlat,zlon), platF, plonF, rd_found_km=rFoundKM, j_prv=0, i_prv=0 )

        #lOK = (jnF>0 and inF>0) ; # a nearest point was found!
        if jnF<0 or inF<0:
            ialive[jP] = 0
            if iverbose>1: print('        ===> I CANCEL this buoy!!! (NO nearest F-point for ',zlat,zlon,')')
                            
        if ialive[jP] == 1:            
            # Ok a nearest point was found!    
            if iverbose>1:
                print('     ==> nearest F-point for ',zlat,zlon,' on NEMO grid:', jnF, inF, '==> lat,lon:',
                      round(platF[jnF,inF],3), round(plonF[jnF,inF],3))
    
            #      o--o           x--o            o--x            o--o
            # 1 => |  | NE   2 => |  | SE    3 => |  | SW    4 => |  | NW
            #      x--o           o--o            o--o            o--x
            iq = Iquadran( (zlat,zlon), platF, plonF, jnF, inF, k_ew_per=-1, lforceHD=True )
            if not iq in [1,2,3,4]:                
                ialive[jP] = 0
                if iverbose>1: print('        ===> I CANCEL this buoy!!! (iQuadran Fuck Up)')
                    
        if ialive[jP] == 1:            
            # Indices of the T-point in the center of the mesh ():
            if   iq==1:
                zjinT[jP,:] = [jnF+1,inF+1]
            elif iq==2:
                zjinT[jP,:] = [jnF  ,inF+1]
            elif iq==3:
                zjinT[jP,:] = [jnF  ,inF]
            elif iq==4:
                zjinT[jP,:] = [jnF+1,inF]
            #
            [jT,iT] = zjinT[jP,:]

            if iverbose>1:                
                print('     =>> coord of center of mesh (T-point) => lat,lon =',
                      round(platT[jT,iT],3), round(plonT[jT,iT],3), '(iq =',iq,')')

            # Test on land-sea mask:
            zmt = maskT[jT,iT] + maskT[jT,iT+1]+maskT[jT+1,iT]+maskT[jT,iT-1]+maskT[jT-1,iT-1]
            if zmt < 5:
                ialive[jP] = 0
                if iverbose>1: print('        ===> I CANCEL this buoy!!! (too close to or over the land-sea mask)')
            
            # Test on sea-ice concentration:
            if len(np.shape(xIceConc))==2:
                #zic = xIceConc[jT,iT]
                zic = 0.2*(xIceConc[jT,iT] + xIceConc[jT,iT+1]+xIceConc[jT+1,iT]+xIceConc[jT,iT-1]+xIceConc[jT-1,iT-1])
                if iverbose>1:
                    #print('     =>> sea-ice concentration =',zic)
                    print('     =>> 5P sea-ice concentration =',zic)
            if zic < rmin_conc:
                ialive[jP] = 0
                if iverbose>1: print('        ===> I CANCEL this buoy!!! (mean 5P ice concentration at T-point of the cell =',zic,')')

        if ialive[jP] == 1:
            # Find the `j,i` indices of the 4 points composing the source mesh that includes the target point
            #  starting with the nearest point
            zJI4vrtc = IDSourceMesh( (zlat,zlon), platF, plonF, jnF, inF, iquadran=iq, k_ew_per=-1, lforceHD=True )
            [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = zJI4vrtc
            # => indexing is anti-clockwize with (j1,i1) beeing the F-point nearest to the buoy...
    
            # Identify the 4 corners (aka F-points) as bottom left, bottom right, upper right & and upper left
            if   iq==1:
                zJIvrtcs[jP,:,:] = [ [ j1, j2, j3, j4 ], [ i1, i2, i3, i4 ] ]
            elif iq==2:
                zJIvrtcs[jP,:,:] = [ [ j2, j3, j4, j1 ], [ i2, i3, i4, i1 ] ]
            elif iq==3:
                zJIvrtcs[jP,:,:] = [ [ j3, j4, j1, j2 ], [ i3, i4, i1, i2 ] ]
            elif iq==4:
                zJIvrtcs[jP,:,:] = [ [ j4, j1, j2, j3 ], [ i4, i1, i2, i3 ] ]
                        
    ### for jP in range(nP)

    return  zjinT, zJIvrtcs










def IsInsideCell( px, py, CellPolygon ):
    '''
         * px, py     : global (not local to the cell) x,y position [m] or [km]
         * CellPolygon: shapely polygon object of a quadrangle mesh constructed with same unit as px and py
    '''
    pnt  = Point(py,px)
    return CellPolygon.contains(pnt)


def ccw( pcA, pcB, pcC ):
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
    return ( ccw(pcA,pcC,pcD) != ccw(pcB,pcC,pcD) ) and ( ccw(pcA,pcB,pcC) != ccw(pcA,pcB,pcD) )


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



