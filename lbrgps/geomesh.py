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
rQarea_max = 800000. ; # max area allowed for Quadrangle [km^2]

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
    ''' The 3 angles of a triangle given the coordinates of the 3 points '''
    from math import sqrt, acos, pi
    #
    nv = len(pCoorT[:,0]) ; # number of vertices...
    if nv != 3:
        print('ERROR [AnglesOfTriangle]: I am only designed for triangles...', nv); exit(0)

    va = [ __distAB2__(pCoorT[i,:], pCoorT[(i+1)%3,:]) for i in range(nv) ]; # Square of the length of each side

    rsa = np.sqrt(va) ; # Length of each side

    # From Cosine law
    vabc = np.array([ acos((va[i]+va[(i+2)%3]-va[(i+1)%3]) / (2*rsa[(i+2)%3]*rsa[i])) for i in range(nv) ])

    return vabc * 180./pi


def AreaOfTriangle(pCoorT):
    ''' Area of a triangle given the coordinates of the 3 points '''
    #
    nv = len(pCoorT[:,0]) ; # number of vertices...
    if nv != 3:
        print('ERROR [AreaOfTriangle]: I am only designed for triangles...') ; exit(0)

    rA = 0.5*(   (pCoorT[0,0]*(pCoorT[1,1] - pCoorT[2,1]))
               + (pCoorT[1,0]*(pCoorT[2,1] - pCoorT[0,1]))
               + (pCoorT[2,0]*(pCoorT[0,1] - pCoorT[1,1])) )
    return rA

def lTisOK( pAngles ):
    '''
    ###  => returns Boolean: True if the triangle has a "decent" shape!
    ###
    ###  * pAngles: the 3 angles of the triangle in degrees 
    '''
    return ( not( np.any(  pAngles   >rTang_max) or np.any( pAngles < rTang_min) ) )


def lQisOK( pAngles, ratio, pArea=None ):
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
    if pArea:
        lOK = (lOK) and (pArea<=rQarea_max)
    if lOK:
        rscore = 1. - np.sum( np.abs(pAngles[:] - 90.) ) / (4.*90.) ; # => 0 if perfect
    return lOK, rscore


def QSpecs2Tri( pTRIAs, it1, it2, pCoor, pnam=[] ):
    '''
    ###  Quadrangle specs from 2 triangles (it1 and it2)
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

    if not it2 in pTRIAs.neighbors[it1,:]:
        print('ERROR: [QSpecs2Tri()] => triangle #'+str(it2)+' is not neighbor with triangle #'+str(it1)+'!'); exit(0)

    [vp1,vp2] = pTRIAs.TriPointIDs[[it1,it2],:] ; # 3 point IDs forming the vertices of triangle #1 and #2

    v2com     = np.intersect1d( vp1, vp2 )
    IDsingle2 = np.setdiff1d(vp2, vp1) ; # Return the ID of the only point of triangle #2 that does not belong to triangle #1
    v2sol     = np.concatenate([np.setdiff1d(vp1, vp2),IDsingle2]) ; # the 2 points  not part of the "common" segment! (1st one belongs to triangle #1)

    jid = IDsingle2[0]
    if ldebug:
        print('   [QSpecs2Tri()] 4 angles of quad when merge triangles #'+str(it1)+' and #'+str(it2)+':')
        print('   [QSpecs2Tri()] ==> the 2 vertices in common: ', v2com, '=', [ pnam[i] for i in v2com ])
        print('   [QSpecs2Tri()] ==> the 2 solitary vertices : ', v2sol, '=', [ pnam[i] for i in v2sol ])        
        print('   [QSpecs2Tri()] ==> point to add to triangle '+str(it1)+' to form a quadrangle is #'+str(jid)+' aka "'+pnam[jid]+'"')

    # Ratio between apparent height and width of the quadrangle
    zcoor_com = np.array([ pCoor[i,:] for i in v2com ])
    zcoor_sol = np.array([ pCoor[i,:] for i in v2sol ])
    Lc11 = sqrt( __distAB2__(zcoor_com[0,:], zcoor_sol[0,:]) ) ; #lilo
    Lc21 = sqrt( __distAB2__(zcoor_com[1,:], zcoor_sol[0,:]) ) ; #lilo
    Lc12 = sqrt( __distAB2__(zcoor_com[1,:], zcoor_sol[1,:]) ) ; #lilo
    Lc22 = sqrt( __distAB2__(zcoor_com[0,:], zcoor_sol[1,:]) ) ; #lilo
    ratio = (Lc11+Lc12) / (Lc21+Lc22)
    del zcoor_com, zcoor_sol, Lc11, Lc21, Lc12, Lc22

    [va1, va2] = pTRIAs.angles()[[it1,it2],:] ; # The 3 angles of each of the 2 triangles:

    # 2 angles associated to the 2 "solitary points":
    va_sol = np.concatenate([ va1[np.where(vp1==v2sol[0])],va2[np.where(vp2==v2sol[1])] ])

    # 2 angles associated for the 2 common points:
    ra1 = np.concatenate([ va1[np.where(vp1==v2com[0])] , va2[np.where(vp2==v2com[0])] ])
    ra2 = np.concatenate([ va1[np.where(vp1==v2com[1])] , va2[np.where(vp2==v2com[1])] ])

    # The four angles of the quadrangle (without any order):
    vIquad, vAquad = np.concatenate([ v2com,v2sol ]), np.concatenate([ [ np.sum(ra1),np.sum(ra2)],va_sol ])

    # Sorting in a counter-clockwize fasion:
    zcoor = np.array( [ pCoor[i,:] for i in vIquad ] )
    vsidx = SortIndicesCCW(zcoor)
    vIquad, vAquad = vIquad[vsidx], vAquad[vsidx]

    if ldebug:
        print('   [QSpecs2Tri()] ==> the 4 angles of the CCW-sorted quadrangle =',vAquad)
        print('   [QSpecs2Tri()]       ===> for ',[ pnam[i] for i in vIquad ])

    # Return the CCW-sorted points and angles for the 4 vertices of the quadrangle:
    return vIquad, vAquad, ratio, pTRIAs.area()[it1]+pTRIAs.area()[it2]


def Tri2Quad( pTRIAs, pCoor, pnam,  iverbose=0 ):
    '''
    ### Attempt to merge triangles into quadrangles:
    ###  Each triangle inside the domain has 3 neighbors, so there are 3 options to merge
    ###  (provided the 3 neighbors are not already merged with someone!)
    '''
    NbT     = pTRIAs.nT
    ivb     = iverbose

    NbQ      = 0
    idxTdead = []   ; # IDs of canceled triangles
    idxTused = []   ; # IDs of triangles already in use in a quadrangle
    Quads    = []

    # Loop along triangles:
    for jT in range(NbT):

        v3pnts  = pTRIAs.TriPointIDs[jT,:] ; # the 3 point IDs composing triangle # jT
        vangles = pTRIAs.angles()[jT]      ; # the 3 angles...
        
        if ivb>0:
            print('\n **************************************************************')
            print(' *** Focus on triangle #'+str(jT)+' =>',[ pnam[i] for i in v3pnts ],'***')
            print(' **************************************************************')

        if not lTisOK(vangles):
            if ivb>0: print('       => disregarding this triangle!!! (an angle >'+str(rTang_max)+' or <'+str(rTang_min)+' degrees!)')
            idxTdead.append(jT) ; # Cancel this triangle

        elif jT in idxTused:
            if ivb>0: print('       => this triangle is in use in an already defined Quad !')


        else:
            #
            # Triangle `jT` has a "decent" shape and has not been used to build a quad yet!
            # -----------------------------------------------------------------------------
            #
            if ivb>1: print('  ==> its 3 angles are:',vangles[:])
            #
            vtmp   = pTRIAs.neighbors[jT,:]
            vnghbs = vtmp[vtmp >= 0] ; # shrink it, only retain non `-1`-flagged values...
            NbN    = len(vnghbs)     ; # number of neighbors
            if ivb>0: print('       => its '+str(NbN)+' neighbor triangles are:', vnghbs)

            NgbrTvalid = [] ; # ID the valid neighbor triangles, i.e.: not dead, not already in use, and decent shape!
            for jN in vnghbs:
                lTok = (not jN in idxTdead)and(not jN in idxTused)and(lTisOK(vangles))
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
                    vidx, vang, rat, area  = QSpecs2Tri( pTRIAs, jT, jN, pCoor) #, pnam=pnam )
                    lQok, score = lQisOK( vang[:], rat, pArea=area ) ;#lilo
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
    del v3pnts, vangles
    
    zQpoints = np.array(Quads)
    zQcoor   = []
    del Quads

    if NbQ>0:

        # Coordinates of the points:
        zQcoor = np.array([ [ pCoor[i,:] for i in zQpoints[jQ,:] ] for jQ in range(NbQ) ]) ; # Shape is (NbQ,4,2) !!!

        # Some sanity checks:
        if len(idxTused)/2 != NbQ or len(idxTused)%2 !=0:
            print('ERROR [Tri2Quad]: agreement between number of merged triangles and created quads!')
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

