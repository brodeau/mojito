#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
##################################################################

from sys import exit
#from os import path, mkdir
import numpy as nmp


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
    iz = nmp.array([0,1,2,3])
    zt, it = OrderCW(xcoor)
    #
    zz[1:4,:] = zt[:0:-1,:]
    iz[1:4]   = it[:0:-1]
    del zt, it
    return zz, iz





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





def QuadAnglesFrom2Tri( pTrgl, pNghb, it1, it2, pCoor, pnam=[] ):
    '''
    ###
    ###          There are `Np` points that define `Nt` triangles !
    ### Input:
    ###
    ###         * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###         * pNghb:    the 3 triangle IDs being the neighbors,   shape: (Nt,3) | origin: `scipy.Delaunay().neighbors`
    ###         * it1, it2: IDs of the two triangles to merge into a quadrangle
    ###         * pCoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    ###         * pnam:   OPTIONAL (DEBUG)  name of each point (string)  , shape: (Np)
    ###
    '''
    ldebug = ( len(pnam)>0 )

    vp0  = pTrgl[it1,:]   ; # 3 point indices forming the triangle

    if not it2 in pNghb[it1,:]:
        print('   [QuadAnglesFrom2Tri()]: ERROR: triangle #'+str(it2)+' is not neighbor with triangle #'+str(it1)+'!'); exit(0)

    vpn  = pTrgl[it2,:]
    v2com = nmp.intersect1d( vp0, vpn )
    v2sol = nmp.concatenate([nmp.setdiff1d(vp0, vpn),nmp.setdiff1d(vpn, vp0)]) ; # First value being the one of triangle it1
    if ldebug:
        print('   [QuadAnglesFrom2Tri()] 4 angles of quad when merge triangles #'+str(it1)+' and #'+str(it2)+':')
        print('   [QuadAnglesFrom2Tri()] ==> the 2 vertices in common: ', v2com, '=', [ pnam[i] for i in v2com ])
        print('   [QuadAnglesFrom2Tri()] ==> the 2 solitary vertices : ', v2sol, '=', [ pnam[i] for i in v2sol ])

    vID_unique_it2 = nmp.setdiff1d(vpn, vp0) ; # Return the unique values in `vpn` that are not in `vp0`.

    jid = vID_unique_it2[0]
    if ldebug:
        print('   [QuadAnglesFrom2Tri()] ==> point to add to triangle '+str(it1)+' to form a quadrangle is #'+str(jid)+' aka "'+pnam[jid]+'"')

    # Now we look at the angles of the 2 triangles: lilo
    va1 = AnglesOfTriangle(it1, pTrgl, pCoor)
    va2 = AnglesOfTriangle(it2, pTrgl, pCoor)
    if ldebug:
        vp1 = pTrgl[it1,:]
        vp2 = pTrgl[it2,:]
        #print('    --- angles for triangle #'+str(it1)+':', va1,'(',[ pnam[i] for i in vp1 ],')')
        #print('    --- angles for triangle #'+str(it2)+':', va2,'(',[ pnam[i] for i in vp2 ],')')

    # 2 angles associated for the 2 solitary points:
    va_sol = nmp.concatenate([ va1[nmp.where(vp0==v2sol[0])],va2[nmp.where(vpn==v2sol[1])] ])
    #if ldebug: print('      --- angles for solo points:', va_sol,'(',[ pnam[i] for i in v2sol ],')')

    # 2 angles associated for the 2 common points:
    ra1 = nmp.concatenate([ va1[nmp.where(vp0==v2com[0])] , va2[nmp.where(vpn==v2com[0])] ])
    #if ldebug: print('      --- angles for common point '+str(v2com[0])+', respect. seen from tri. #'+str(it1)+' & #'+str(it2)+':',ra1)
    ra2 = nmp.concatenate([ va1[nmp.where(vp0==v2com[1])] , va2[nmp.where(vpn==v2com[1])] ])
    #if ldebug: print('      --- angles for common point '+str(v2com[1])+', respect. seen from tri. #'+str(it1)+' & #'+str(it2)+':',ra2)

    # The four angles of the quadrangle (without any order):
    vIquad = nmp.concatenate([     v2com     ,     v2sol              ])
    vAquad = nmp.concatenate([ [ nmp.sum(ra1), nmp.sum(ra2) ], va_sol ])

    # Sorting in a counter-clockwize fasion:
    zcoor = nmp.array( [ pCoor[i,:] for i in vIquad ] )
    vsidx = SortIndicesCCW(zcoor)
    vIquad = vIquad[vsidx]
    vAquad = vAquad[vsidx]

    if ldebug:
        print('   [QuadAnglesFrom2Tri()] ==> the 4 angles of the CCW-sorted quadrangle =',vAquad)
        print('   [QuadAnglesFrom2Tri()]       ===> for ',[ pnam[i] for i in vIquad ])

    # Return the CCW-sorted points and angles for the 4 vertices of the quadrangle:
    return vIquad, vAquad


def lQuadOK( pAngles ):
    '''
    ###     Tells if a uqdrangle should be kept or rejected (outrageously unrectangular shape)
    ###       + gives a score [0-1] (for now based on
    ###     TO DO: add the height/width ratio !!!! (we want somethibg closer to as square than a rectangle)
    ###
    ### Input:
    ###
    ###         * pAngles: 1D array containing the 4 vertex angles of the quadrangle [degrees]
    ###
    '''
    rscore = -1.
    lOK = ( not( any(pAngles>140.) or any(pAngles<30.) ) )
    if lOK:
        rscore = 1. - nmp.sum( nmp.abs(pAngles[:] - 90.) ) / (4.*90.) ; # => 0 if perfect

    return lOK, rscore



def Triangles2Quadrangle( pTrgl, pNghb, it1, it2, pCoor, pnam=[] ):
    '''
    ###
    ###          There are `Np` points that define `Nt` triangles !
    ### Input:
    ###
    ###         * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###         * pNghb:    the 3 triangle IDs being the neighbors,   shape: (Nt,3) | origin: `scipy.Delaunay().neighbors`
    ###         * it1, it2: IDs of the two triangles to merge into a quadrangle
    ###         * pCoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
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
        print('                                    ===>' ,[ (round(pCoor[i,1],2),round(pCoor[i,0],2)) for i in quad ] )


    Zcoor = nmp.array([[pCoor[i,0],pCoor[i,1]] for i in quad ])

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

def ThreeAngles(pCoor):
    from math import sqrt, acos, pi
    ## pCoor: `x,y` coordinates of the 3 points shape=(3,2)
    # Square of lengths be a2, b2, c2
    a2 = __lengthSquare__(pCoor[1,:], pCoor[2,:])
    b2 = __lengthSquare__(pCoor[0,:], pCoor[2,:])
    c2 = __lengthSquare__(pCoor[0,:], pCoor[1,:])

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



def AnglesOfTriangle(kT, pTrgl, pCoor):
    '''
    #### Wrapper for AnglesOfTriangle
    ###
    ###  *   kT :     ID of triangle we are working with
    ###  * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###  * pCoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    '''
    v3p = pTrgl[kT,:] ; # 3 point IDs composing the triangle with ID `kt`
    zcoorT = nmp.array([ pCoor[j,:] for j in v3p ])
    return ThreeAngles( zcoorT )

def lTriangleOK(kT, pTrgl, pCoor):
    '''
    #### Wrapper for AnglesOfTriangle
    ###  => returns Boolean
    ###
    ###  *   kT :     ID of triangle we are working with
    ###  * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###  * pCoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    '''
    v3p = pTrgl[kT,:] ; # 3 point IDs composing the triangle with ID `kt`
    zcoorT = nmp.array([ pCoor[j,:] for j in v3p ])
    va = ThreeAngles( zcoorT )
    return (nmp.any(va>110.) or nmp.any(va<30.))



####################################

def Triangles2Quads( pTrgl, pNghb, pCoor, pnam,  iverbose=0 ):
    '''
    ### Attempt to merge triangles into quadrangles:
    ###  Each triangle inside the domain has 3 neighbors, so there are 3 options to merge
    ###  (provided the 3 neighbors are not already merged with someone!)
    '''
    (NbT,_) = nmp.shape(pTrgl)
    ivb     = iverbose

    NbQ      = 0
    idxTdead = []   ; # IDs of canceled triangles
    idxTused = []   ; # IDs of triangles already in use in a quadrangle
    Quads    = []

    # Loop along triangles:
    for jT in range(NbT):

        v3pnts = pTrgl[jT,:] ; # the 3 point IDs composing triangle # jT

        if ivb>0:
            print('\n **************************************************************')
            print(' *** Focus on triangle #'+str(jT)+' =>',[ pnam[i] for i in v3pnts ],'***')
            print(' **************************************************************')

        if lTriangleOK(jT, pTrgl, pCoor):
            if ivb>0: print('       => disregarding this triangle!!! (an angle >120. or <30 degrees!)')
            idxTdead.append(jT) ; # Cancel this triangle

        elif jT in idxTused:
            if ivb>0: print('       => this triangle is in use in an already defined Quad !')


        else:
            # Triangle `jT` has a "decent" shape and has not been used to build a quad yet!
            # -----------------------------------------------------------------------------
            #
            vtmp   = pNghb[jT,:]
            vnghbs = vtmp[vtmp >= 0] ; # shrink it, only retain non `-1`-flagged values...
            NbN    = len(vnghbs)     ; # number of neighbors
            if ivb>0: print('       => its '+str(NbN)+' neighbor triangles are:', vnghbs)

            NgbrTvalid = [] ; # ID the valid neighbor triangles, i.e.: not dead, not already in use, and decent shape!
            for jN in vnghbs:
                lTok = (not jN in idxTdead)and(not jN in idxTused)and(not lTriangleOK(jN, pTrgl, pCoor))
                if lTok:
                    NgbrTvalid.append(jN)
                    if ivb>1: print('          ==> triangle '+str(jN)+' is valid!')
                else:
                    if ivb>1: print('          ==> triangle '+str(jN)+' is NOT valid!')

            if len(NgbrTvalid)>0:
                # `jT` is a valid+available triangle with at least one valid neighbor in `NgbrTvalid`
                #   => need to check which of the neighbors in `NgbrTvalid` gives the best quadrangle!
                if ivb>0: print('       => valid neighbors for triangle #'+str(jT)+':',NgbrTvalid)

                #NNTok  = 0  ; # number of neighbor triangles are ok for forming a "decent" quad
                vjNok  = [] ; # stores the neighbor triangles that are ok for forming a "decent" quad
                vscore = [] ; #   => store their score of "OK-ness" !
                xidx   = []
                for jN in NgbrTvalid:
                    if ivb>1: print('          ==> trying neighbor triangle '+str(jN)+':')
                    vidx, vang  = QuadAnglesFrom2Tri( pTrgl, pNghb, jT, jN, pCoor) #, pnam=pnam )
                    lQok, score = lQuadOK( vang[:] )
                    cc = 'does NOT'
                    if lQok:
                        vjNok.append(jN)
                        vscore.append(score)
                        xidx.append(vidx)
                        cc = 'does'
                    if ivb>1: print('            ===> "triangles '+str(jT)+'+'+str(jN)+'" '+cc+' give a valid Quad!')
                # Now we have to chose the best neighbor triangle to use (based on the score):
                if len(vjNok)>0:
                    if ivb>0: print('       => We have '+str(len(vjNok))+' Quad candidates!')
                    xidx = nmp.array(xidx)
                    iwin = nmp.argmax(vscore)
                    jN   = vjNok[iwin] ; # our winner neighbor triangle
                    Quads.append(xidx[iwin,:]) ; # saving the winner Quad!
                    idxTused.append(jT)
                    idxTused.append(jN)
                    NbQ = NbQ+1
                    if ivb>0: print('         ==> Selected Quad: "triangles '+str(jT)+'+'+str(jN)+'" give Quad',[pnam[i] for i in xidx[iwin,:]])

            else:
                if ivb>0: print('       => No valid neighbors for this triangle...')
    ## -- for jT in range(NbT) --
    zQuads = nmp.array(Quads)
    del Quads

    # Some sanity checks:
    if len(idxTused)/2 != NbQ or len(idxTused)%2 !=0:
        print('ERROR [Triangles2Quads]: agreement between number of merged triangles and created quads!')
        exit(0)

    if ivb>0: print('\n *** SUMMARY ***')
    if ivb>1:
        print('       => Triangles sucessfully merged into "acceptable" Quads:')
        print('       ==>',idxTused,)
    if ivb>0:
        print('       => Summary about the '+str(NbQ)+' Quads generated:')
        for jQ in range(NbQ):
            print('        * Quad #'+str(jQ)+' => ', zQuads[jQ,:], '(', [ pnam[i] for i in zQuads[jQ,:] ],')')
        print('')

    return zQuads


