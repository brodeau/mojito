#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import exit
#from os import path, mkdir
import numpy as np

# Sloppy:
rTang_min =  15. ; # minimum angle tolerable in a triangle [degree]
rTang_max = 115. ; # maximum angle tolerable in a triangle [degree]
rQang_min =  65. ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 120. ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.5 ; # value that `1 - abs(L/H)` should not overshoot!

# Strict:
#rTang_min =  20. ; # minimum angle tolerable in a triangle [degree]
#rTang_max = 100. ; # maximum angle tolerable in a triangle [degree]
#rQang_min =  75. ; # minimum angle tolerable in a quadrangle [degree]
#rQang_max = 105. ; # maximum angle tolerable in a quadrangle [degree]
#rdRatio_max = 0.25 ; # value that `1 - abs(L/H)` should not overshoot!


def __distAB2__(pC1, pC2):
    ''' Square of the distance between 2 points based on their [x,y] coordinates '''
    xDiff = pC1[0] - pC2[0]
    yDiff = pC1[1] - pC2[1]
    return xDiff*xDiff + yDiff*yDiff


def AnglesOfTriangle(pCoorT):
    from math import sqrt, acos, pi
    ## pCoorT: `x,y` coordinates of the 3 points shape=(3,2)

    nv = len(pCoorT[:,0]) ; # number of vertices...
    if not nv in [3,4]:
        print('ERROR [AnglesOfTriangle]: I am designed for triangles or quadrangles...')
        exit(0)

    va = [ __distAB2__(pCoorT[i,:], pCoorT[(i+1)%3,:]) for i in range(nv) ]; # Square of the length of each side

    rsa = np.sqrt(va) ; # Length of each side

    # From Cosine law
    vabc = np.array([ acos((va[i]+va[(i+2)%3]-va[(i+1)%3]) / (2*rsa[(i+2)%3]*rsa[i])) for i in range(nv) ])

    return vabc * 180./pi



def QuadSpecsFrom2Tri( pTrgl, pNghb, it1, it2, pCoor, pnam=[] ):
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
    from math  import sqrt
    from .util import SortIndicesCCW

    ldebug = ( len(pnam)>0 )

    if not it2 in pNghb[it1,:]:
        print('ERROR: [QuadSpecsFrom2Tri()] => triangle #'+str(it2)+' is not neighbor with triangle #'+str(it1)+'!')
        exit(0)

    vp0 = pTrgl[it1,:] ; # 3 point IDs forming the vertices of triangle #1
    vpn = pTrgl[it2,:] ; # 3 point IDs forming the vertices of triangle #2

    v2com     = np.intersect1d( vp0, vpn )
    IDsingle2 = np.setdiff1d(vpn, vp0) ; # Return the ID of the only point of triangle #2 that does not belong to triangle #1
    v2sol     = np.concatenate([np.setdiff1d(vp0, vpn),IDsingle2]) ; # the 2 points that are not part of the "common" segment! (1st value belongs tp triangle #1)
    if ldebug:
        print('   [QuadSpecsFrom2Tri()] 4 angles of quad when merge triangles #'+str(it1)+' and #'+str(it2)+':')
        print('   [QuadSpecsFrom2Tri()] ==> the 2 vertices in common: ', v2com, '=', [ pnam[i] for i in v2com ])
        print('   [QuadSpecsFrom2Tri()] ==> the 2 solitary vertices : ', v2sol, '=', [ pnam[i] for i in v2sol ])


    jid = IDsingle2[0]
    if ldebug:
        print('   [QuadSpecsFrom2Tri()] ==> point to add to triangle '+str(it1)+' to form a quadrangle is #'+str(jid)+' aka "'+pnam[jid]+'"')

    # Ratio between apparent height and width of the quadrangle
    zcoor_com = np.array([ pCoor[i,:] for i in v2com ])
    zcoor_sol = np.array([ pCoor[i,:] for i in v2sol ])
    Lc11 = sqrt( __distAB2__(zcoor_com[0,:], zcoor_sol[0,:]) ) ; #lilo
    Lc21 = sqrt( __distAB2__(zcoor_com[1,:], zcoor_sol[0,:]) ) ; #lilo
    Lc12 = sqrt( __distAB2__(zcoor_com[1,:], zcoor_sol[1,:]) ) ; #lilo
    Lc22 = sqrt( __distAB2__(zcoor_com[0,:], zcoor_sol[1,:]) ) ; #lilo
    ratio = (Lc11+Lc12) / (Lc21+Lc22)
    del zcoor_com, zcoor_sol, Lc11, Lc21, Lc12, Lc22

    # Now we look at the angles of the 2 triangles: lilo
    va1 = __triangle_angles__( pCoor, vp0 )
    va2 = __triangle_angles__( pCoor, vpn )
    if ldebug:
        vp1 = pTrgl[it1,:]
        vp2 = pTrgl[it2,:]
        #print('    --- angles for triangle #'+str(it1)+':', va1,'(',[ pnam[i] for i in vp1 ],')')
        #print('    --- angles for triangle #'+str(it2)+':', va2,'(',[ pnam[i] for i in vp2 ],')')

    # 2 angles associated for the 2 solitary points:
    va_sol = np.concatenate([ va1[np.where(vp0==v2sol[0])],va2[np.where(vpn==v2sol[1])] ])
    #if ldebug: print('      --- angles for solo points:', va_sol,'(',[ pnam[i] for i in v2sol ],')')

    # 2 angles associated for the 2 common points:
    ra1 = np.concatenate([ va1[np.where(vp0==v2com[0])] , va2[np.where(vpn==v2com[0])] ])
    #if ldebug: print('      --- angles for common point '+str(v2com[0])+', respect. seen from tri. #'+str(it1)+' & #'+str(it2)+':',ra1)
    ra2 = np.concatenate([ va1[np.where(vp0==v2com[1])] , va2[np.where(vpn==v2com[1])] ])
    #if ldebug: print('      --- angles for common point '+str(v2com[1])+', respect. seen from tri. #'+str(it1)+' & #'+str(it2)+':',ra2)

    # The four angles of the quadrangle (without any order):
    vIquad = np.concatenate([     v2com     ,     v2sol              ])
    vAquad = np.concatenate([ [ np.sum(ra1), np.sum(ra2) ], va_sol ])

    # Sorting in a counter-clockwize fasion:
    zcoor = np.array( [ pCoor[i,:] for i in vIquad ] )
    vsidx = SortIndicesCCW(zcoor)
    vIquad = vIquad[vsidx]
    vAquad = vAquad[vsidx]

    if ldebug:
        print('   [QuadSpecsFrom2Tri()] ==> the 4 angles of the CCW-sorted quadrangle =',vAquad)
        print('   [QuadSpecsFrom2Tri()]       ===> for ',[ pnam[i] for i in vIquad ])

    # Return the CCW-sorted points and angles for the 4 vertices of the quadrangle:
    return vIquad, vAquad, ratio


def lQuadOK( pAngles, ratio ):
    '''
    ###     Tells if a quadrangle should be kept or rejected (outrageously unrectangular shape)
    ###       + gives a score [0-1] (for now based on
    ###     TO DO: add the height/width ratio !!!! (we want somethibg closer to as square than a rectangle)
    ###
    ### Input:
    ###         * pAngles: 1D array containing the 4 vertex angles of the quadrangle [degrees]
    ###         * ratio:   ratio between apparent length and height of the quadrangle (a square would give 1)
    ###
    '''
    rscore = -1.
    lOK = ( not( np.any(pAngles>rQang_max) or np.any(pAngles<rQang_min) or abs(1.-ratio)>rdRatio_max ) )
    if lOK:
        rscore = 1. - np.sum( np.abs(pAngles[:] - 90.) ) / (4.*90.) ; # => 0 if perfect

    return lOK, rscore



def __triangle_angles__( pCoor, p3p ):
    '''
    ###  *   kT :     ID of triangle we are working with
    ###  * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###  * pCoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    '''
    zcoorT = np.array([ pCoor[j,:] for j in p3p ])
    return AnglesOfTriangle( zcoorT )

def lTriangleOK(kT, pTrgl, pCoor):
    '''
    ###  => returns Boolean: True if the triangle has a "decent" shape!
    ###
    ###  *   kT :     ID of triangle we are working with
    ###  * pTrgl:    the 3 point IDs    forming the triangles, shape: (Nt,3) | origin: `scipy.Delaunay().simplices`
    ###  * pCoor:    (lon,lat) coordinates of each point,      shape: (Np,2)
    '''
    v3p = pTrgl[kT,:] ; # 3 point IDs composing the triangle with ID `kt`
    zcoorT = np.array([ pCoor[j,:] for j in v3p ])
    va = AnglesOfTriangle( zcoorT )
    lOK = ( not( np.any(  va   >rTang_max) or np.any( va < rTang_min) ) )
    #
    return lOK



####################################

def Triangles2Quads( pTrgl, pNghb, pCoor, pnam,  iverbose=0 ):
    '''
    ### Attempt to merge triangles into quadrangles:
    ###  Each triangle inside the domain has 3 neighbors, so there are 3 options to merge
    ###  (provided the 3 neighbors are not already merged with someone!)
    '''
    (NbT,_) = np.shape(pTrgl)
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

        if not lTriangleOK(jT, pTrgl, pCoor):
            if ivb>0: print('       => disregarding this triangle!!! (an angle >'+str(rTang_max)+' or <'+str(rTang_min)+' degrees!)')
            idxTdead.append(jT) ; # Cancel this triangle

        elif jT in idxTused:
            if ivb>0: print('       => this triangle is in use in an already defined Quad !')


        else:
            # DEBUG: interested what are the angles:
            if ivb>1: print('  ==> its 3 angles are:',__triangle_angles__(pCoor, v3pnts))
            # DEBUG.

            # Triangle `jT` has a "decent" shape and has not been used to build a quad yet!
            # -----------------------------------------------------------------------------
            #
            vtmp   = pNghb[jT,:]
            vnghbs = vtmp[vtmp >= 0] ; # shrink it, only retain non `-1`-flagged values...
            NbN    = len(vnghbs)     ; # number of neighbors
            if ivb>0: print('       => its '+str(NbN)+' neighbor triangles are:', vnghbs)

            NgbrTvalid = [] ; # ID the valid neighbor triangles, i.e.: not dead, not already in use, and decent shape!
            for jN in vnghbs:
                lTok = (not jN in idxTdead)and(not jN in idxTused)and(lTriangleOK(jN, pTrgl, pCoor))
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
                    vidx, vang, rat  = QuadSpecsFrom2Tri( pTrgl, pNghb, jT, jN, pCoor) #, pnam=pnam )
                    lQok, score = lQuadOK( vang[:], rat )
                    cc = 'does NOT'
                    if lQok:
                        vjNok.append(jN)
                        vscore.append(score)
                        xidx.append(vidx)
                        cc = 'does'
                    if ivb>1:
                        print('            ===> "triangles '+str(jT)+'+'+str(jN)+'" '+cc+' give a valid Quad!')
                        if not lQok: print('              ====> the 4 angles + ratio:',vang[:],rat)

                # Now we have to chose the best neighbor triangle to use (based on the score):
                if len(vjNok)>0:
                    if ivb>0: print('       => We have '+str(len(vjNok))+' Quad candidates!')
                    xidx = np.array(xidx)
                    iwin = np.argmax(vscore)
                    jN   = vjNok[iwin] ; # our winner neighbor triangle
                    Quads.append(xidx[iwin,:]) ; # saving the winner Quad!
                    idxTused.append(jT)
                    idxTused.append(jN)
                    NbQ = NbQ+1
                    if ivb>0: print('         ==> Selected Quad: "triangles '+str(jT)+'+'+str(jN)+'" give Quad',[pnam[i] for i in xidx[iwin,:]])

            else:
                if ivb>0: print('       => No valid neighbors for this triangle...')
    ## -- for jT in range(NbT) --

    zQpoints = np.array(Quads)
    zQcoor   = []
    del Quads

    if NbQ>0:

        # Coordinates of the points:
        zQcoor = np.array([ [ pCoor[i,:] for i in zQpoints[jQ,:] ] for jQ in range(NbQ) ]) ; # Shape is (NbQ,4,2) !!!

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
                print('        * Quad #'+str(jQ)+' => ', zQpoints[jQ,:], '(', [ pnam[i] for i in zQpoints[jQ,:] ],')')

    else:
        print('\n WARNING => No Quads could be generated! :(')
    print('')

    return zQpoints, zQcoor


