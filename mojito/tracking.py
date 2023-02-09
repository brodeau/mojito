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


def CrossedEdge( pP1, pP2, j4vert, i4vert, pY, pX,  iverbose=0):
    '''
        * pP1    : point 1 as [y1,x1]
        * pP2    : point 2 as [y2,x2]
        * j4vert : vector of length 4, containg the 4 j-indices of the 4 vertices of the mesh
        * i4vert : vector of length 4, containg the 4 i-indices of the 4 vertices of the mesh
    '''
    [ jbl, jbr, jur, jul ] = j4vert[:]
    [ ibl, ibr, iur, iul ] = i4vert[:]
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

