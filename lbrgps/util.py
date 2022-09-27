#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import exit
#from os import path, mkdir
import numpy as nmp
#from re import split
#import climporn as cp

#idebug = 0


def LoadDist2CoastNC( cNCfile ):
    from climporn import chck4f
    from netCDF4  import Dataset
    #
    print('\n *** [util.LoadDist2CoastNC()] Loading "distance to coast" from file:')
    print('      '+cNCfile)
    chck4f(cNCfile)
    with Dataset(cNCfile) as id_in:
        vlon  = id_in.variables['lon'][:]
        vlat  = id_in.variables['lat'][:]
        xdist = id_in.variables['dist'][:,:]
        print('       => ok!\n')
        #
    return vlon, vlat, xdist


def Dist2Coast( lon0, lat0, plon, plat, pdist2coat ):
    '''
       Returns the distance to the nearest coast of a given point (lon,lat)
        INPUT:
          * lon0, lat0: coordinates (scalars) [degrees East], [degrees West]
          * plon:       1D array of longitudes assosiated to `pdist2coat`
          * plat:       1D array of latitudes  assosiated to `pdist2coat`
          * pdist2coat: 2D array containing rasterized distance to coast (originally read in NC file) [km]
    '''
    from climporn import degE_to_degWE
    #
    rx = nmp.mod( lon0, 360. ) ; # [0:360] frame
    vx = nmp.mod( plon, 360. ) ; # [0:360] frame
    # Are we dangerously close to the [0--360] cut?
    #  => then will work in the [-180:180] frame:
    if rx > 355.:
        rx = degE_to_degWE( rx )
        vx = degE_to_degWE( vx )
        #print(' lon0, lat0 =', rx, lat0)
    ip = nmp.argmin( nmp.abs(  vx[:] - rx  ) )
    jp = nmp.argmin( nmp.abs(plat[:] - lat0) )
    #print(' ip, jp =', ip, jp)
    #print(' Nearest lon, lat =', plon[ip], plat[jp])
    del vx, rx
    return max( pdist2coat[jp,ip] , 0. )



def OrderCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''

    # sort the points based on their x-coordinates
    isortX  = nmp.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]
    print('LOLO: Sorted by longitude => isortX =', isortX)

    
    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]
    
    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = nmp.argsort(leftMost[:,1])
    leftMost   = leftMost[isortL,:]
    (tl, bl)   = leftMost
    [i1l, i2l] = isortML[isortL]
    print('LOLO: isortL, idxL =', isortL, [i1l, i2l] )

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = nmp.argsort(rghtMost[:,1])
    rghtMost   = rghtMost[isortR,:]
    (tr, br)   = rghtMost
    [i1r, i2r] = isortMR[isortR]
    print('LOLO: isortR, idxR =',isortR, [i1r, i2r] )
    
    #idx = nmp.concatenate([idxL,idxR])
    #isort = nmp.array([isortX[i] for i in nmp.concatenate([isortL,isortR])])
    #print('LOLO, idx=', idx)

    del isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR
    
    # return the coordinates in top-left, top-right,
    # bottom-right, and bottom-left order
    return nmp.array([tl, tr, br, bl], dtype="float32") , nmp.array([i1l, i1r, i2r, i2l])


def OrderCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''
    zz = xcoor.copy()
    zz[1:4,:] = OrderCW(xcoor)[:0:-1,:]
    return zz, 





def SortIndicesCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter-clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''

    # sort the points based on their x-coordinates
    isortX  = nmp.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]
    #print('LOLO: Sorted by longitude => isortX =', isortX)

    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]
    
    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = nmp.argsort(leftMost[:,1])
    [i1l, i2l] = isortML[isortL]
    #print('LOLO: isortL, idxL =', isortL, [i1l, i2l] )

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = nmp.argsort(rghtMost[:,1])
    [i1r, i2r] = isortMR[isortR]
    #print('LOLO: isortR, idxR =',isortR, [i1r, i2r] )
    
    #idx = nmp.concatenate([idxL,idxR])
    #isort = nmp.array([isortX[i] for i in nmp.concatenate([isortL,isortR])])
    #print('LOLO, idx=', idx)

    del xSorted, isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR
    
    return nmp.array([i1l, i1r, i2r, i2l])





def WouldBeValidQuad( pTrgl, pNghb, it1, it2, pcoor, pnam=[] ):
    '''
    ###
    ###          There are `Np` points that define `Nt` triangles !
    ### Input:
    ###        
    ###         * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###         * pNghb:    the 3 triangle IDs being the neighbors,   shape: (Nt,3) | origin: `scipy.Delaunay().neighbors`
    ###         * it1, it2: IDs of the two triangles to merge into a quadrangle
    ###         * pcoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    ###         * pnam:   OPTIONAL (DEBUG)  name of each point (string)  , shape: (Np)
    ###
    '''
    ldebug = ( len(pnam)>0 )
    
    vp0  = pTrgl[it1,:]   ; # 3 point indices forming the triangle

    if not it2 in pNghb[it1,:]:
        print('   [util.WouldBeValidQuad()]: ERROR: triangle #'+str(it2)+' is not neighbor with triangle #'+str(it1)+'!'); exit(0)

    vpn  = pTrgl[it2,:]
    v2com = nmp.intersect1d( vp0, vpn )
    v2sol = nmp.concatenate([nmp.setdiff1d(vp0, vpn),nmp.setdiff1d(vpn, vp0)]) ; # First value being the one of triangle it1
    if ldebug:
        print('   [util.WouldBeValidQuad()] Triangles #'+str(it1)+' and #'+str(it2)+':')
        print('   [util.WouldBeValidQuad()] ==> the 2 vertices in common: ', v2com, '=', [ pnam[i] for i in v2com ])
        print('   [util.WouldBeValidQuad()] ==> the 2 solitary vertices : ', v2sol, '=', [ pnam[i] for i in v2sol ])

    vID_unique_it2 = nmp.setdiff1d(vpn, vp0) ; # Return the unique values in `vpn` that are not in `vp0`.

    jid = vID_unique_it2[0]
    if ldebug:
        print('   [util.WouldBeValidQuad()] ==> point to add to triangle '+str(it1)+' to form a quadrangle is #'+str(jid)+' aka "'+pnam[jid]+'"')

    # Now we look at the angles of the 2 triangles: lilo
    va1 = AnglesOfTriangle(it1, pTrgl, pcoor)
    va2 = AnglesOfTriangle(it2, pTrgl, pcoor)
    if ldebug:
        vp1 = pTrgl[it1,:]
        vp2 = pTrgl[it2,:]
        print('    --- angles for triangle #'+str(it1)+':', va1,'(',[ pnam[i] for i in vp1 ],')')
        print('    --- angles for triangle #'+str(it2)+':', va2,'(',[ pnam[i] for i in vp2 ],')')

    # 2 angles associated to v2sol:
    idsolo = v2sol[0]
    ([ii],) = nmp.where(vp0==idsolo)
    print(' angle for solo point '+str(idsolo)+' ('+pnam[idsolo]+') =',va1[ii])
    
    idsolo = v2sol[1]
    ([ii],) = nmp.where(vpn==idsolo)
    print(' angle for solo point '+str(idsolo)+' ('+pnam[idsolo]+') =',va2[ii])



    
    for iv in v2com:
        print('   - lolo: vertex#'+str(iv)+' ('+pnam[iv]+'):')




    
    exit(0)
        
    quad = nmp.concatenate( [ vID_unique_it2, vp0 ] ) ; # This is our quadrangle !!!

    iOrder = [0,1,2,3]
    if ldebug:
        print('   [util.WouldBeValidQuad()] ==> order after triangle merge:', iOrder, '=>',[ pnam[i]         for i in quad ] )
        print('                                    ===>' ,[ (round(pcoor[i,1],2),round(pcoor[i,0],2)) for i in quad ] )


    Zcoor = nmp.array([[pcoor[i,0],pcoor[i,1]] for i in quad ])

    # Ordering the 4 points in a clockwise fashion:
    iOrder = SortIndicesCCW(Zcoor)

    # Reordering quad:
    quadCCW = quad[iOrder]

    del quad, Zcoor
    
    return quadCCW





def Triangles2Quadrangle( pTrgl, pNghb, it1, it2, pcoor, pnam=[] ):
    '''
    ###
    ###          There are `Np` points that define `Nt` triangles !
    ### Input:
    ###        
    ###         * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###         * pNghb:    the 3 triangle IDs being the neighbors,   shape: (Nt,3) | origin: `scipy.Delaunay().neighbors`
    ###         * it1, it2: IDs of the two triangles to merge into a quadrangle
    ###         * pcoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    ###         * pnam:   OPTIONAL (DEBUG)  name of each point (string)  , shape: (Np)
    ###
    '''
    ldebug = ( len(pnam)>0 )
    
    vp0  = pTrgl[it1,:]   ; # 3 point indices forming the triangle

    if not it2 in pNghb[it1,:]: print('   [util.Triangles2Quadrangle()]: ERROR: triangle #'+str(it2)+' is not neighbor with triangle #'+str(it1)+'!'); exit(0)

    vpn  = pTrgl[it2,:]
    v2com = nmp.intersect1d( vp0, vpn )
    if ldebug:
        print('   [util.Triangles2Quadrangle()] ==> the 2 vertices in common between triangles #'+str(it1)+' and #'+str(it2)+': ', v2com, '=', [ pnam[i] for i in v2com ])

    vID_unique_it2 = nmp.setdiff1d(vpn, vp0) ; # Return the unique values in `vpn` that are not in `vp0`.

    jid = vID_unique_it2[0]
    if ldebug:
        print('   [util.Triangles2Quadrangle()] ==> point to add to triangle '+str(it1)+' to form a quadrangle is #'+str(jid)+' aka "'+pnam[jid]+'"')

    quad = nmp.concatenate( [ vID_unique_it2, vp0 ] ) ; # This is our quadrangle !!!

    iOrder = [0,1,2,3]
    if ldebug:
        print('   [util.Triangles2Quadrangle()] ==> order after triangle merge:', iOrder, '=>',[ pnam[i]         for i in quad ] )
        print('                                    ===>' ,[ (round(pcoor[i,1],2),round(pcoor[i,0],2)) for i in quad ] )


    Zcoor = nmp.array([[pcoor[i,0],pcoor[i,1]] for i in quad ])

    # Ordering the 4 points in a clockwise fashion:
    iOrder = SortIndicesCCW(Zcoor)

    # Reordering quad:
    quadCCW = quad[iOrder]

    del quad, Zcoor
    
    return quadCCW





def __lengthSquare__(X, Y):
    xDiff = X[0] - Y[0]
    yDiff = X[1] - Y[1]
    return xDiff*xDiff + yDiff*yDiff

def ThreeAngles(pcoor):
    from math import sqrt, acos, pi
    ## pcoor: `x,y` coordinates of the 3 points shape=(3,2)
    # Square of lengths be a2, b2, c2
    a2 = __lengthSquare__(pcoor[1,:], pcoor[2,:])
    b2 = __lengthSquare__(pcoor[0,:], pcoor[2,:])
    c2 = __lengthSquare__(pcoor[0,:], pcoor[1,:])
    
    # length of sides be a, b, c
    a = sqrt(a2);
    b = sqrt(b2);
    c = sqrt(c2);
    
    # From Cosine law
    alpha = acos((b2 + c2 - a2) / (2 * b * c));
    betta = acos((a2 + c2 - b2) / (2 * a * c));
    gamma = acos((a2 + b2 - c2) / (2 * a * b));
    
    # Converting to degree
    return nmp.array([ alpha, betta, gamma ])* 180. / pi
    


def AnglesOfTriangle(kT, pTrgl, pcoor):
    '''
    #### Wrapper for AnglesOfTriangle
    ###
    ###  *   kT :     ID of triangle we are working with
    ###  * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###  * pcoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    '''    
    v3p = pTrgl[kT,:] ; # 3 point IDs composing the triangle with ID `kt`    
    zcoorT = nmp.array([ pcoor[j,:] for j in v3p ])    
    return ThreeAngles( zcoorT )
    
def AnglesOfTriangleNotOK(kT, pTrgl, pcoor):
    '''
    #### Wrapper for AnglesOfTriangle
    ###  => returns Boolean
    ###
    ###  *   kT :     ID of triangle we are working with
    ###  * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###  * pcoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    '''    
    v3p = pTrgl[kT,:] ; # 3 point IDs composing the triangle with ID `kt`    
    zcoorT = nmp.array([ pcoor[j,:] for j in v3p ])    
    va = ThreeAngles( zcoorT )
    return (nmp.any(va>110.) or nmp.any(va<30.))
    
