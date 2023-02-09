#import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#from math import atan2,pi


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


def CrossedEdge( pP1, pP2, j4vert, i4vert, pY, pX,  iverbose=0 ):
    '''
        * pP1    : point 1 as [y1,x1]
        * pP2    : point 2 as [y2,x2]
        * j4vert : vector of length 4, containg the 4 j-indices of the 4 vertices of the mesh
        * i4vert : vector of length 4, containg the 4 i-indices of the 4 vertices of the mesh
    '''
    [ ibl, ibr, iur, iul ] = i4vert[:]
    [ jbl, jbr, jur, jul ] = j4vert[:]
    #
    for kk in range(4):
        kp1 = (kk+1)%4
        j1 = j4vert[kk]
        i1 = i4vert[kk]
        j2 = j4vert[kp1]
        i2 = i4vert[kp1]
        #
        ll = intersect2Seg( pP1, pP2,  [pY[j1,i1],pX[j1,i1]], [pY[j2,i2],pX[j2,i2]] )
        if ll: break
        #
    if iverbose>0:
        vdir = ['bottom', 'right-hand', 'upper', 'left-hand']
        print('    [CrossedEdge()]: particle is crossing the '+vdir[kk]+' edge of the mesh!')
    return kk+1




def NewHostCell( kcross, pP1, pP2, j4vert, i4vert, pY, pX,  iverbose=0 ):
    '''
       Find in which adjacent cell, the point has moved into !
       
       * kcross : integer describing which (old) cell edge the point has crossed

       In rare cases, the buoy could possibly move into diagonally adjacent cells
       (not only bottom,right,upper,left cells...)
       These "diagonally-adjacent" meshes yields output: 5,6,7,8
    '''
    #
    [ jbl, jbr, jur, jul ] = j4vert[:]
    [ ibl, ibr, iur, iul ] = i4vert[:]
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



def UpdtInd4NewCell( knhc, j4vert, i4vert, kjT, kiT ):
    
    if   knhc==1:
        # went down:
        j4vert[:] = j4vert[:] - 1
        kjT = kjT-1
    elif knhc==5:
        print('LOLO: WE HAVE A 5 !!!!')
        # went left + down
        i4vert[:] = i4vert[:] - 1
        j4vert[:] = j4vert[:] - 1
        kiT = kiT-1
        kjT = kjT-1
    elif knhc==6:
        print('LOLO: WE HAVE A 6 !!!!')
        # went right + down
        i4vert[:] = i4vert[:] + 1
        j4vert[:] = j4vert[:] - 1
        kiT = kiT+1
        kjT = kjT-1                    
    elif knhc==2:
        # went to the right
        i4vert[:] = i4vert[:] + 1
        kiT = kiT+1
    elif knhc==3:
        # went up:
        j4vert[:] = j4vert[:] + 1
        kjT = kjT+1
    elif knhc==7:
        print('LOLO: WE HAVE A 7 !!!!')
        # went right + up
        i4vert[:] = i4vert[:] + 1
        j4vert[:] = j4vert[:] + 1
        kiT = kiT+1
        kjT = kjT+1
    elif knhc==8:
        print('LOLO: WE HAVE A 8 !!!!')
        # went right + up
        i4vert[:] = i4vert[:] - 1
        j4vert[:] = j4vert[:] + 1
        kiT = kiT-1
        kjT = kjT+1
    elif knhc==4:
        # went to the left
        i4vert[:] = i4vert[:] - 1
        kiT = kiT-1
    else:
        print('ERROR: unknown direction, knhc=',knhc)
        exit(0)

    return j4vert, i4vert, kjT, kiT

