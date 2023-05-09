#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import exit
from os import path
import numpy as np


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



def __distAB2__(pC1, pC2):
    ''' Square of the distance between 2 points based on their [x,y] coordinates '''
    xDiff = pC1[0] - pC2[0]
    yDiff = pC1[1] - pC2[1]
    return xDiff*xDiff + yDiff*yDiff


def LengthsOfTriangle(pCoorT):
    ''' The length of the 3 segments defining a triangle given the coordinates of the 3 points '''
    from math import sqrt, acos, pi
    #
    nv = len(pCoorT[:,0]) ; # number of vertices...
    if nv != 3:
        print('ERROR [LengthsOfTriangle]: I am only designed for triangles...', nv); exit(0)

    va = [ __distAB2__(pCoorT[i,:], pCoorT[(i+1)%3,:]) for i in range(nv) ]; # Square of the length of each side
    return np.sqrt(va)


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


def AreaOfQuadrangle( pX4, pY4 ):
    '''
    === Area of a quadrangle from the (x,y) coordinates of its 4 vertices ===
           (see for example Bouchat et al. 2022, Eq.5 )
    # Input:
             * pX4: the 4 x-coordinates (as a length, not longitude!)
             * pY4: the 4 y-coordinates (as a length, not latitude!)
    # Output:
             * the area of the quadrangle
    '''
    nv = len(pX4[:]) ; # number of vertices...
    if len(pX4[:]) != 4 or len(pY4[:]) != 4:
        print('ERROR [AreaOfQuadrangle]: I am only designed for quadrangles! => 4 vertices!!!') ; exit(0)
    #
    rA = np.sum( np.array([ pX4[k]*pY4[(k+1)%4] - pX4[(k+1)%4]*pY4[k] for k in range(4) ]) )
    #
    return 0.5*rA

def QuadsAreas( pQVcoor ):
    '''
    # Input:
             * pQVcoor: [shape = (nQ,4,2)] the (x,y) coordinates of the 4 vertices of the `nQ` quads
    # Output:
             * the areas of the quadrangles
    '''
    shp = np.shape(pQVcoor)
    if len(shp)!=3:
        print('ERROR [QuadsAreas]: wrong number of dims for `pQVcoor`! =>',len(shp),' we need 3!'); exit(0)
    (nq,n4,n2) = shp
    if n4!=4 or n2!=2:
        print('ERROR [QuadsAreas]: wrong shape for `pQVcoor`!'); exit(0)
    zA = np.zeros(nq) - 9999.
    for i in range(nq):
        zA[i] = AreaOfQuadrangle( pQVcoor[i,:,0] , pQVcoor[i,:,1])
    #
    return zA




def lTisOK( pAngles, pArea=None, anglR=(15.,115.), areaR=(0.,5.e20) ):
    '''
    ###  => returns Boolean: True if the triangle has a "decent" shape!
    ###
    ###  * pAngles: the 3 angles of the triangle in degrees
    '''
    lOK = ( not( np.any(pAngles>anglR[1]) or np.any(pAngles<anglR[0]) ) )
    if lOK and pArea:
        lOK = ( (pArea>areaR[0]) and (pArea<=areaR[1]) )
    return lOK


def lQisOK( pAngles, ratio, pArea=None, ratioD=0.5, anglR=(65.,120.), areaR=(0.,8.e20), areaIdeal=None ):
    '''
    ###     Tells if a quadrangle should be accepted or rejected (outrageously unrectangular shape)
    ###       + gives a score [0-1] (for now based on
    ###     TO DO: add the height/width ratio !!!! (we want somethibg closer to as square than a rectangle)
    ###
    ### Input:
    ###         * pAngles: 1D array containing the 4 vertex angles of the quadrangle [degrees]
    ###         * ratio: max ratio between apparent length and height of the quadrangle
    ###                  =>  >= 1 !!! (a square would give 1)
    ###
    '''
    #print('LOLO [lQisOK]: A, A_min, A_max, A_ideal =',int(pArea), int(areaR[0]), int(areaR[1]), int(areaIdeal), end=' ')
    rscore = -1.
    lOK = ( pArea>0. )
    if lOK:
        lOK = ( not( np.any(pAngles>anglR[1]) or np.any(pAngles<anglR[0]) or ratio-1.>ratioD ) )
    if lOK and pArea:
        lOK = ( pArea>areaR[0] and pArea<areaR[1] )
    if lOK:
        if areaIdeal:
            # We base the score on the deviation from the ideal area:
            rscore = 1. / (1. + 0.001*np.abs(pArea - areaIdeal))
            #print('==> rscore (area) =',round(rscore,4))
        else:
            # We base the score on the angles of the quadrangle:
            rscore = 1. - np.sum( np.abs(pAngles[:] - 90.) ) / (4.*90.) ; # => 0 if perfect
            #print('==> rscore (angles) =',round(rscore,4))
    #print('')
    return lOK, rscore



def QSpecsFrom2T( p3Pnt, pAngl, it1, it2, pCoor, pnam=[] ):
    '''
    ###  Quadrangle specs obtained from "would-merge" the 2 triangles with IDs `it1` and `it2`
    ###
    ###       Note: about the shape pf pCoor and pnam:
    ###             a cloud of `nP` points defines `nT` triangles ! => #fixme? (is it bad to pass the entire `pCoor` ?)
    ###
    ### Input:
    ###         * p3Pnt:   the 3 point IDs forming the triangles, shape: (2,3) | origin: `scipy.Delaunay().simplices`
    ###         * pAngl:   the 3 angles of the triangle,          shape: (2,3)
    ###         * it1, it2: IDs of the two triangles to merge into a quadrangle
    ###         * pCoor:  coordinates of each point,      shape: (nP,2)
    ###         * pnam:   OPTIONAL (DEBUG)  name of each point (string)  , shape: (nP)
    ###
    '''
    from math  import sqrt
    from .util import SortIndicesCCW

    ldebug = ( len(pnam)>0 )

    [vp1, vp2] = p3Pnt[:,:] ; # 3 point IDs forming the vertices of triangle #1 and #2
    [va1, va2] = pAngl[:,:] ; # The 3 angles of each of the 2 triangles:

    v2com     = np.intersect1d( vp1, vp2 )
    IDsingle2 = np.setdiff1d(vp2, vp1) ; # Return the ID of the only point of triangle #2 that does not belong to triangle #1
    v2sol     = np.concatenate([np.setdiff1d(vp1, vp2),IDsingle2]) ; # the 2 points  not part of the "common" segment!
    #                                                                   => 1st one belongs to triangle #1

    jid = IDsingle2[0]
    if ldebug:
        print('   [QSpecsFrom2T()] 4 angles of quad when merge triangles #'+str(it1)+' and #'+str(it2)+':')
        print('   [QSpecsFrom2T()] ==> the 2 vertices in common: ', v2com, '=', [ pnam[i] for i in v2com ])
        print('   [QSpecsFrom2T()] ==> the 2 solitary vertices : ', v2sol, '=', [ pnam[i] for i in v2sol ])
        print('   [QSpecsFrom2T()] ==> point to add to triangle '+str(it1)+' to form a quadrangle is #'+str(jid)+' aka "'+pnam[jid]+'"')

    # Ratio between apparent height and width of the quadrangle
    zcoor_com = np.array([ pCoor[i,:] for i in v2com ])
    zcoor_sol = np.array([ pCoor[i,:] for i in v2sol ])
    Lc11 = sqrt( __distAB2__(zcoor_com[0,:], zcoor_sol[0,:]) )
    Lc21 = sqrt( __distAB2__(zcoor_com[1,:], zcoor_sol[0,:]) )
    Lc12 = sqrt( __distAB2__(zcoor_com[1,:], zcoor_sol[1,:]) )
    Lc22 = sqrt( __distAB2__(zcoor_com[0,:], zcoor_sol[1,:]) )
    rh1 = Lc11+Lc12
    rh2 = Lc21+Lc22
    ratio = max( rh1/rh2 , rh2/rh1 ) ; # it must be >= 1. !!! (square =>1, >1 rectangle-like regardless of the direction)
    del zcoor_com, zcoor_sol, Lc11, Lc21, Lc12, Lc22, rh1, rh2

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
    zcoor[:] = zcoor[vsidx]

    # Area of the would-be quadrangle:
    zarea = AreaOfQuadrangle( zcoor[:,0], zcoor[:,1] )
    
    if ldebug:
        print('   [QSpecsFrom2T()] ==> area of the would-be quadrangle =',zarea)
        print('   [QSpecsFrom2T()] ==> the 4 angles of the CCW-sorted quadrangle =',vAquad)
        print('   [QSpecsFrom2T()]       ===> for ',[ pnam[i] for i in vIquad ])

    # Return the CCW-sorted points indices and angles for the 4 vertices of the quadrangle, the area of the quads, as well as ratio
    return vIquad, vAquad, zarea, ratio


def TriPntIDs2QuaPntIDs( xPntID ):
    '''
        * xPntID: [ (nQ,4) array of integers] the 4 point IDs composing the Quad,
                  point IDs are those based on the nT points defining the initial triangles
                  from which Quads were built (nQ << nT)

     => RETURNS:
        *
        * zPntID: same as xPntID but point IDs are those based on only the nQ points defining the Quads
    '''
    zPntIDs1D = np.unique( xPntID.flatten() )
    nP = len(zPntIDs1D) ; # number of points that defines the nQ quadrangles
    zPntID = xPntID.copy()
    zPntID[:,:] = -1
    for jP in range(nP):
        ii = zPntIDs1D[jP]
        idx = np.where(xPntID==ii)
        zPntID[idx] = jP

    return np.array( zPntID, dtype=int )


def Tri2Quad( pTRIAs, anglRtri=(15.,115.), ratioD=0.5, anglR=(65.,120.), areaR=(0.,8.e20), areaIdeal=None, idbglev=0 ):
    '''
    ### Attempt to merge triangles into quadrangles:
    ###  Each triangle inside the domain has 3 neighbors, so there are 3 options to merge
    ###  (provided the 3 neighbors are not already merged with someone!)
    ###
    ### RETURNS:
    ## zPxy : [nP,2] array of floats] the coordinates of the nP points that define the Quads
    ## zPIDs : [nP] vector of int] IDs of the points of zPxy (IDs left from the original cloud of points on which triangles are build)
    ## zPtime : [nP] vector of floats] the time of the points of zPxy
    ## zQPQ   : [nQ,4] array of int] the 4 point indices composing the quad, in counter-clockwize
    ## zQname: [nQ]   vector of strings]  a string to identify each quadrangle
    '''
    NbT = pTRIAs.nT
    ivb = idbglev

    NbQ      = 0
    idxTdead = []   ; # IDs of canceled triangles
    idxTused = []   ; # IDs of triangles already in use in a quadrangle
    Quads    = []   ; # Valid quads identified...

    zXY  = pTRIAs.PointXY.copy()      ; # shape: (nP,2)
    zIDs   = pTRIAs.PointIDs.copy()     ; # shape: (nP)
    ztime  = pTRIAs.PointTime.copy()     ; # shape: (nP)
    zcN    = pTRIAs.PointNames.copy()   ; # shape: (nP)
    Z3Pnts = pTRIAs.MeshVrtcPntIdx.copy() ; # shape: (nT,3)
    Znghbs = pTRIAs.NeighborIDs.copy()  ; # shape: (nT,3)
    Zangls = pTRIAs.angles().copy()     ; # shape: (nT,3)
    Zareas = pTRIAs.area().copy()       ; # shape: (nT)

    # Loop along triangles:
    for jT in range(NbT):

        v3pnts  = Z3Pnts[jT,:] ; # the 3 point IDs composing triangle # jT
        vangles = Zangls[jT,:] ; # the 3 angles...
        rarea   = Zareas[jT]

        if ivb>1:
            print('\n **************************************************************')
            print(' *** Focus on triangle #'+str(jT)+' =>',[ zcN[i] for i in v3pnts ],'***')
            print(' **************************************************************')

        # What is an acceptable area range for a triangle (to become part of a  areaR[0] <quadrangle <areaR[1])
        #   => let's be tolerant: `> areaR[0]/3` & `< areaR[1]/1.4` 
        if not lTisOK(vangles, pArea=rarea, anglR=anglRtri, areaR=(areaR[0]/3.,areaR[1]/1.4)):
            if ivb>1: print('       => disregarding this triangle!!! (because of extreme angles or area)')
            idxTdead.append(jT) ; # Cancel this triangle

        elif jT in idxTused:
            if ivb>1: print('       => this triangle is in use in an already defined Quad !')

        else:
            # Triangle `jT` has a "decent" shape and has not been used to build a quad yet
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if ivb>1: print('  ==> its 3 angles are:',vangles[:])

            vtmp   = Znghbs[jT,:]
            vnghbs = vtmp[vtmp >= 0] ; # shrink it, only retain non `-1`-flagged values...
            if ivb>1: print('       => its '+str(len(vnghbs))+' neighbor triangles are:', vnghbs)

            NgbrTvalid = [] ; # ID the valid neighbor triangles, i.e.: not dead, not already in use [REMOVED:, and decent shape]!
            for jN in vnghbs:
                lTok = ( (not jN in idxTdead) and (not jN in idxTused) )
                if lTok:
                    NgbrTvalid.append(jN)
                    if ivb>1: print('          ==> triangle '+str(jN)+' is valid!')
                else:
                    if ivb>1: print('          ==> triangle '+str(jN)+' is NOT valid!')

            if len(NgbrTvalid)>0:
                # `jT` is a valid+available triangle with at least one valid neighbor in `NgbrTvalid`
                #   => need to check which of the neighbors in `NgbrTvalid` gives the best quadrangle!
                if ivb>1: print('       => valid neighbors for triangle #'+str(jT)+':',NgbrTvalid)

                vjNok  = [] ; # stores the neighbor triangles that are ok for forming a "decent" quad
                vscore = [] ; #   => store their score of "OK-ness" !
                xPids  = []
                for jN in NgbrTvalid:
                    if ivb>1: print('          ==> trying neighbor triangle '+str(jN)+':')

                    vInd4V, vAng4V, zA, rat  = QSpecsFrom2T( Z3Pnts[[jT,jN],:], Zangls[[jT,jN],:], jT, jN, zXY, pnam=[] )
                    
                    lQok, score = lQisOK( vAng4V, rat, pArea=zA, ratioD=ratioD, anglR=anglR, areaR=areaR, areaIdeal=areaIdeal )

                    cc = 'does NOT'
                    if lQok:
                        vjNok.append(jN)
                        vscore.append(score)
                        xPids.append(vInd4V)
                        cc = 'does'
                    if ivb>1:
                        print('            ===> "triangles '+str(jT)+'+'+str(jN)+'" '+cc+' give a valid Quad!')
                        if not lQok: print('              ====> the 4 angles + ratio:',vAng4V[:],rat)

                # Now we have to chose the best neighbor triangle to use (based on the score):
                if len(vjNok)>0:
                    if ivb>1: print('       => We have '+str(len(vjNok))+' Quad candidates!')
                    xPids = np.array(xPids)
                    iwin = np.argmax(vscore)
                    jN   = vjNok[iwin] ; # our winner neighbor triangle
                    Quads.append(xPids[iwin,:]) ; # saving the winner Quad!
                    idxTused.append(jT)
                    idxTused.append(jN)
                    NbQ = NbQ+1
                    if ivb>1: print('   ==> New Quad: "triangles '+str(jT)+'+'+str(jN)+'" => ',[zcN[i] for i in xPids[iwin,:]])

            else:
                if ivb>1: print('       => No valid neighbors for this triangle...')

        if ivb>1: print('**************************************************************')
    ## -- for jT in range(NbT) --
    del Z3Pnts, Zangls, Zareas, Znghbs, v3pnts, vangles

    zQPT = np.array(Quads)
    del Quads

    # Following arrays might be returned void if no valid Quad is identified...
    zPxy  = [] ; 
    zPIDs   = [] ; # 
    zPtime  = [] ; # 
    zQpind  = [] ; # the 4 point indices composing the quad, counter-clockwize 
    zQname = [] ; # 

    if NbQ>0:
        # Coordinates of the points in use by Quads!
        zvPntIdx = np.unique( zQPT.flatten() ) ; # isolates the point IDs that are in use by the identified+valid Quads...
        zPxy = np.array([ zXY[i,:] for i in zvPntIdx ])
        zPIDs  = np.array([  zIDs[i]   for i in zvPntIdx ], dtype=int )
        zPtime = np.array([ ztime[i]   for i in zvPntIdx ], dtype=int )
        del zvPntIdx


        zQname = np.array( [ zcN[zQPT[jQ,0]]+'-'+zcN[zQPT[jQ,1]]+'-'+zcN[zQPT[jQ,2]]+'-'+zcN[zQPT[jQ,3]] for jQ in range(NbQ) ],
                            dtype='U32' )

        # Point IDs (from original triangle cloud) are now translated to the points that remains for Quads:
        zQpind = TriPntIDs2QuaPntIDs(zQPT)

        # Some sanity checks:
        if len(idxTused)/2 != NbQ or len(idxTused)%2 !=0:
            print('ERROR [Tri2Quad]: agreement between number of merged triangles and created quads!')
            exit(0)

        if ivb>1: print('\n *** SUMMARY ***')
        if ivb>1:
            print('       => Triangles sucessfully merged into "acceptable" Quads:')
            print('       ==>',idxTused,)
        if ivb>1:
            print('       => Summary about the '+str(NbQ)+' Quads generated:')
            for jQ in range(NbQ):
                print('        * Quad #'+str(jQ)+' => ', zQPT[jQ,:], '('+zQname[jQ]+')')
                if ivb>2:
                    vx = np.array([ zPxy[i,0] for i in zQpind[jQ,:] ])
                    vy = np.array([ zPxy[i,1] for i in zQpind[jQ,:] ])
                    print('     X-coor:', vx[:])
                    print('     Y-coor:', vy[:])

    else:
        print('\n WARNING => No Quads could be generated! :(')
    print('')

    # Check if some "crossed" quads:
    zQcoor = np.array( [ zPxy[zQpind[jQ,:],:] for jQ in range(NbQ) ]) ; #magic, the (x,y) coordinates of the 4 vert. of each Quad! shape:(nQ,4,2) !
    for jQ in range(NbQ):
        if intersect2Seg( zQcoor[jQ,0,:], zQcoor[jQ,1,:], zQcoor[jQ,2,:], zQcoor[jQ,3,:] ):
            print('WARNING [Tri2Quad()]: we have crossed fucked-up quadrangle!!! => name =',zQname[jQ])
            # Correction: p3 moved to pos#1 (apparently no other bad combinations arise...):
            zQcoor[jQ,1:,:] = [ zQcoor[jQ,3,:], zQcoor[jQ,1,:], zQcoor[jQ,2,:] ]
            zQpind[jQ,1:]   = [ zQpind[jQ,3],   zQpind[jQ,1],   zQpind[jQ,2]   ]
            zQname[jQ]      = zcN[zQPT[jQ,0]]+'-'+zcN[zQPT[jQ,3]]+'-'+zcN[zQPT[jQ,1]]+'-'+zcN[zQPT[jQ,2]]
            if intersect2Seg( zQcoor[jQ,0,:], zQcoor[jQ,1,:], zQcoor[jQ,2,:], zQcoor[jQ,3,:] ):
                print('ERROR [Tri2Quad()]: fucked-up Quad! Do not know how to fix!'); exit(0)
            print('WARNING [Tri2Quad()]: fixed! => new name is: ',zQname[jQ])

    
    if ivb>0:
        # DEBUG: identify negative areas??? How did they pass through `lQisOK`? #fixme
        zQcoor = np.array( [ zPxy[zQpind[jQ,:],:] for jQ in range(NbQ) ]) ; #magic, the (x,y) coordinates of the 4 vert. of each Quad!
        zareas = QuadsAreas( zQcoor )
        if np.any(zareas<0.):
            print('ERROR [Tri2Quad]: at least a quad with negative area!')
            (idxFU,) = np.where(zareas<0.)
            print('   for Quads with following indices:',idxFU)
            print(' Area(s) => ', zareas[idxFU])
            print(' Name(s) => ', zQname[idxFU])
            print('   => fix `mjt.Tri2Quad` so these do not go through!!!')
            exit(0)
        
    return zPxy, zPIDs, zPtime, zQpind, zQname





def QuadVrtcTDev( pPcoor, pPids, pTime, pQpnts, pQnam, t_dev_cancel, mode='stdev' ):
    '''
    #
    # Get rid of quadrangles that have excessively asynchronous vertices:
    #  * t_dev_cancel: maximum time deviation allowed accros the 4 vertices of a quadrangle [s]
    #
    '''
    from .util import StdDev
    #
    if not mode in ['stdev','maxdiff']:
        print('ERROR [QuadVrtcTDev()]: mode "'+mode+'" is unknown!'); exit(0)
    
    (nQ,_) = np.shape(pQpnts)
    zQVtime = np.array([ pTime[i] for i in pQpnts ]) ; # => same shape as `pQpnts` !
    idxKeep = []
    std_max = 0.
    for jQ in range(nQ):
        z4t = zQVtime[jQ,:]
        if mode=='stdev':
            zstd = StdDev( np.mean(z4t), z4t )
            if zstd < t_dev_cancel:
                idxKeep.append(jQ)
            #if zstd>std_max:
            #    std_max=zstd
        #print(' * [QuadVrtcTDev()]: max point-time-position StDev found within a Quadrangles is =',round(std_max/60.,2),'min')
        #
        elif mode=='maxdiff':
            diff_max = np.max(z4t) - np.min(z4t)
            if diff_max < t_dev_cancel:
                idxKeep.append(jQ)

    idxKeep = np.array( idxKeep , dtype=int )
    nQn = len(idxKeep)
    if nQn<=0:
        print(' * [QuadVrtcTDev()]: 0 quads left!!! => exiting!!!'); exit(0)
    if nQn<nQ:
        # Some Quads must be cancelled!!!!
        print(' * [QuadVrtcTDev()]: we shall disregard '+str(nQ-nQn)+' quads /'+str(nQ)+
              ' because excessive time dev. accros the 4 vertices (>'+str(int(t_dev_cancel/60.))+' min)')
        nQ = nQn
        zQpnts = np.zeros((nQ,4),dtype=int)
        zQnam  = np.zeros( nQ,   dtype='U32')
        zQpnts[:,:] = pQpnts[idxKeep,:] ; #=> not enough because it's indices in the current n. of points => must be corrected later!
        zQnam[:]    = pQnam[idxKeep]
        # Now arrays with `nP` as a dimension:
        zPidx = np.unique( zQpnts ); # =>indices !!!
        zPids = pPids[zPidx]
        (nP,) = np.shape(zPids)
        zTime  = np.zeros( nP   , dtype=int)
        zPcoor = np.zeros((nP,2))
        jP = 0
        for jid in zPids:
            zQpnts[np.where(zQpnts==zPidx[jP])] = jP
            ([idx],)     = np.where(pPids==jid)
            zTime[jP]    =  pTime[idx]
            zPcoor[jP,:] = pPcoor[idx,:]
            jP+=1
        return zPcoor, zPids, zTime, zQpnts, zQnam
    else:        
        return pPcoor, pPids, pTime, pQpnts, pQnam ; # Nothing canceled!




def PDVfromPos( pdt, pXY1, pXY2, pA1, pA2,  xtime1=[], xtime2=[], iverbose=0 ):
    '''
        Computes spatial (x,y) partial derivatives of the velocity vector.
        The velocity vector is constructed from the two consecutive X,Y positions
        given as `pXY1, pXY2`...
          --- `nq` is the number of quadrangles provided ---
        * pdt        : time interval between the 2 consecutive position time [s]
        * pXY1, pXY2 : shape=(nq,4,2), the two consecutive X,Y positions of the 4 vertices of each quad [km]
        * pA1,  pA2  : shape=(nq),     the two consecutive areas of each quad [km^2]
        * ( xtime1, xtime2 ): shape=(nq,4), time [s] for each point of each quadrangle

       Note???: keeping calculated velocities in km/s rather is in theory not a problem because, since A is in km^2, 
                in computation of strain rates these 2 cancel each other...

    '''
    #
    (nq,n4,n2) = np.shape(pXY1)
    if n4!=4 or n2!=2: print('ERROR [PDVfromPos()]: wrong shape for `pXY1`!'); exit(0)
    if np.shape(pXY2)!=(nq,n4,n2): print('ERROR [PDVfromPos()]: `pXY1` & `pXY2` do not agree in shape!'); exit(0)
    if np.shape(pA1 )!=(nq,) or np.shape(pA2)!=(nq,): print('ERROR [PDVfromPos()]: wrong shape for `pA1` or `pA2`!'); exit(0)

    
    # Velocities at the 4 vertices at center of time interval:
    if np.shape(xtime1)==(nq,n4) and np.shape(xtime2)==(nq,n4):
        print('   * [PDVfromPos()]: using accurate time of points for `dt` in estimate of velocities!')                        
        zdt = xtime2.copy() * 0.
        zdt = xtime2[:,:]-xtime1[:,:]
        #
        if np.any(zdt==0.):
            from .util import epoch2clock
            (idxBad,_) = np.where(zdt==0.)
            for jq in idxBad:
                print('  ==> jq, zdt[jq,0], xtime1[jq,0], xtime2[jq,0]=', zdt[jq,0], epoch2clock(xtime1[jq,0]), epoch2clock(xtime2[jq,0]) )
            print(' ERROR [PDVfromPos]: `xtime1` and `xtime2` are identical for some points!!! (see above!)')
            exit(0)
        #
        zU = np.array( [ (pXY2[:,k,0] - pXY1[:,k,0])/zdt[:,k] for k in range(4) ] ).T ; # 1000 because X,Y in km !!!
        zV = np.array( [ (pXY2[:,k,1] - pXY1[:,k,1])/zdt[:,k] for k in range(4) ] ).T ; # 1000 because X,Y in km !!!
        del zdt
    else:
        # Using provided `pdt` rather than exact time difference between points:
        zU = np.array( [  pXY2[:,k,0] - pXY1[:,k,0]           for k in range(4) ] ).T / pdt ; # 1000 because X,Y in km !!!
        zV = np.array( [  pXY2[:,k,1] - pXY1[:,k,1]           for k in range(4) ] ).T / pdt ; # 1000 because X,Y in km !!!

    # Positions of the 4 vertices of each quad at center of time interval:
    zX , zY = np.zeros((nq,4)) , np.zeros((nq,4))
    zX = 0.5*( pXY1[:,:,0] + pXY2[:,:,0] )
    zY = 0.5*( pXY1[:,:,1] + pXY2[:,:,1] )

    # Area of quadrangles at center of time interval:
    zA = np.zeros(nq) - 999.
    #zA[:] = 0.5*( pA1 + pA2 ) ; #* 1.e6 ; # 1.e6 => from km^2 to m^2
    zA[:] = pA1 ; #* 1.e6 ; # 1.e6 => from km^2 to m^2
    if iverbose>1:
        print('')
        for jQ in range(0,nq,100):
            print('  areas =',np.round(zA[jQ],3),'km^2, U =',np.round(zU[jQ,:],5),'m/s, V =',np.round(zV[jQ,:],5),'m/s')
    
    # Partial derivatives:
    #  --- the fact that units for coordinates was km and for area km^2 has no importance because it cancels out,
    #      we are looking to something in [s-1]
    zdUdxy = np.zeros((nq,2))
    zdVdxy = np.zeros((nq,2))
    for jQ in range(nq):
        zd = 1./(2.*zA[jQ])
        zdUdxy[jQ,0] =  np.sum( np.array([ (zU[jQ,(k+1)%4] + zU[jQ,k]) * (zY[jQ,(k+1)%4] - zY[jQ,k]) for k in range(4) ]) ) * zd
        zdUdxy[jQ,1] = -np.sum( np.array([ (zU[jQ,(k+1)%4] + zU[jQ,k]) * (zX[jQ,(k+1)%4] - zX[jQ,k]) for k in range(4) ]) ) * zd
        zdVdxy[jQ,0] =  np.sum( np.array([ (zV[jQ,(k+1)%4] + zV[jQ,k]) * (zY[jQ,(k+1)%4] - zY[jQ,k]) for k in range(4) ]) ) * zd
        zdVdxy[jQ,1] = -np.sum( np.array([ (zV[jQ,(k+1)%4] + zV[jQ,k]) * (zX[jQ,(k+1)%4] - zX[jQ,k]) for k in range(4) ]) ) * zd

    if iverbose>1:
        for jQ in range(0,nq,100):
            print('  dU/dx =',zdUdxy[jQ,0],'1/s, dU/dy =',zdUdxy[jQ,1],'1/s')
            print('  dV/dx =',zdVdxy[jQ,0],'1/s, dV/dy =',zdVdxy[jQ,1],'1/s\n')

    return zX, zY, 1000.*zU, 1000.*zV, zdUdxy, zdVdxy, zA


def DivPDV( pdUdxy, pdVdxy ):
    # Divergence aka invariant SigmaI
    return pdUdxy[:,0] + pdVdxy[:,1]

def ShearPDV( pdUdxy, pdVdxy ):
    # Maximum Shear Strain Rate aka invariant SigmaII
    ztp1 = pdUdxy[:,0] - pdVdxy[:,1]
    ztp2 = pdUdxy[:,1] + pdVdxy[:,0]
    return np.sqrt( ztp1*ztp1 + ztp2*ztp2 )

def Shear2PDV( pdUdxy, pdVdxy ):
    # Square of the Maximum Shear Strain Rate aka invariant SigmaII
    ztp1 = pdUdxy[:,0] - pdVdxy[:,1]
    ztp2 = pdUdxy[:,1] + pdVdxy[:,0]
    return ztp1*ztp1 + ztp2*ztp2





def rJIrJJtoCoord( pJJs, pJIs, pIDs, plon_t, plon_u, plat_t, plat_v ):
    '''
        Get `ji,jj` from trackmice trajectories and convert it to Cartesian coordinates in km
    '''
    #
    (nB,) = np.shape(pIDs)
    if len(pJJs)!=nB or len(pJIs)!=nB:
        print('ERROR [rJIrJJtoCoord()]: len(pJJs)!=nB or len(pJIs)!=nB'); exit(0)
    zlon, zlat = np.zeros(nB), np.zeros(nB)
    
    for jb in range(nB):
        rjj, rji = pJJs[jb], pJIs[jb]

        jj, rj = int(rjj)-1, rjj%1.   ; # F2C !
        ji, ri = int(rji)-1, rji%1.   ; # F2C !
        # #fixme: I'm still not sure whether 0.5 means T point or U,V point !!! => find out !!!
        ####      => for now, assume T is at 0 and U is at 0.5 (that might be the opposite...)

        #  --  working with geographic coordinates rather than cartesian coordinates...
        if ri <= 0.5:
            # 0<=ri<=0.5 ==> then we must interpolate between T_i and U_i:
            rlon = 2.*(0.5-ri)*plon_t[jj,ji] + 2.*ri*plon_u[jj,ji]
        else:
            # 0.5<ri<1 ==> then we must interpolate between U_i and T_i+1:
            rlon = 2.*(1.-ri)*plon_u[jj,ji] + 2*(ri-0.5)*plon_t[jj,ji+1]

        if rj <= 0.5:
            # 0<=rj<=0.5 ==> then we must interpolate between T_j and V_j:
            rlat = 2.*(0.5-rj)*plat_t[jj,ji] + 2.*rj*plat_v[jj,ji]
        else:
            # 0.5<rj<1 ==> then we must interpolate between V_j and T_j+1:
            rlat = 2.*(1.-rj)*plat_v[jj,ji] + 2*(rj-0.5)*plat_t[jj+1,ji]

        zlon[jb] = rlon
        zlat[jb] = rlat

    # Conversion from Geo coordinates lon,lat to Cartesian `NorthPolarStereo` projection in [km]
    from cartopy.crs import PlateCarree, NorthPolarStereo
    crs_src = PlateCarree()
    crs_trg = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    zx,zy,_ = crs_trg.transform_points(crs_src, zlon, zlat).T

    zCCoor = np.array( [  zx ,  zy  ] ).T / 1000. ; #  Cartesian coordinates [km]
    zGCoor = np.array( [ zlon, zlat ] ).T         ; # Geographic coordinates [degrees]
    
    return zGCoor, zCCoor







def MaskCoastal( pGC, mask=[], rMinDistLand=100, fNCdist2coast='dist2coast_4deg_North.nc', convArray='F' ):
    '''
        * pGC: shape(n,2) => coordinates in the form [ lon, lat ] if convArray='F' or [ lat, lon ] if convArray='C' !
        * rMinDistLand: minimum distance to coast allowed [km]
        * fNCdist2coast   : netCDF file containing "distance to coast" info
    
        RETURNS: mask==0 for buoys that must be deleted....
    '''
    from .ncio import LoadDist2CoastNC
    from .util import Dist2Coast
    #
    if (not path.exists(fNCdist2coast)):
        print('ERROR [MaskCoastal()]: provide '+fNCdist2coast+' does not exist!!!'); exit(0)        
    if rMinDistLand<=0:
        print('ERROR [MaskCoastal()]: rMinDistLand<=0 !!!'); exit(0)

    (nB,_) = np.shape(pGC)
    zcoor = np.zeros(np.shape(pGC))    
    if convArray == 'C':
        zcoor[:,:] = np.array( [ pGC[:,1], pGC[:,0] ] ).T
    else:
        zcoor[:,:] = np.array( pGC )
    
    mask1d = np.zeros(nB, dtype='i1') + 1
    if np.shape(mask) == (nB,):
        mask1d[:] = mask[:]            
        
    vlon_dist, vlat_dist, xdist = LoadDist2CoastNC( fNCdist2coast ) ; # Load `distance to coast` data...

    for jb in range(nB):
        if mask1d[jb]==1:
            rD = Dist2Coast( zcoor[jb,0], zcoor[jb,1], vlon_dist, vlat_dist, xdist )
            if rD < rMinDistLand:
                mask1d[jb] = 0
    
    nBn = np.sum(mask1d)
    print('   +++ [MaskCoastal()]: found '+str(nB-nBn)+' buoys to remove due to excessive shore proximity!')
    
    return mask1d





def ShrinkArrays( pmask, pNam, pIDs, pGC, pXY, ptime, recAxis=2 ):
    '''
        RETURNS: shrinked version of input arrays
                 => all elements corresponding to points where
                    mask==0 are deleted!
       * pmask, pNam, pIDs => shape = (nP)
       * pGC, pXY          => shape = (nP,2,nrec)

    '''
    cEM = 'ERROR [geomesh.ShrinkArrays()]:'
    
    nBi = len(pmask) ; # number of initial points
    
    if len(pIDs)!=nBi or len(pNam)!=nBi:
        print(cEM+' shape problem => `len(pIDs)!=nBi or len(pNam)!=nBi` !!!,',nBi, len(pIDs),len(pNam)); exit(0)
    if recAxis==2:
        (nP,nd,nrec) = np.shape(pXY)
        (nt1,ntr)    = np.shape(ptime)
    elif recAxis==0:
        (nrec,nP,nd) = np.shape(pXY)
        (ntr,nt1)    = np.shape(ptime)
    if nP!=nBi or nd!=2: 
        print(cEM+' shape problem => `nP!=nBi or nd!=2` !!!, ',nBi,nP,nd); exit(0)
    if np.shape(pXY)!=np.shape(pGC):
        print(cEM+' shape problem => `shape(pXY)!=shape(pGC)` !!!'); exit(0)                        
    if (nt1,ntr)!=(nP,nrec):
        print(cEM+' shape problem => `shape(ptime)!=(nP,nrec)` !!!'); exit(0)
        
    nBo  = np.sum(pmask) ; # number of points to keep
    zIDs = np.zeros( nBo, dtype=int  )
    zNam = np.zeros( nBo, dtype='U32')
    if recAxis==2:
        zGC = np.zeros((nBo,2,nrec))
        zXY = np.zeros((nBo,2,nrec))
        ztm = np.zeros((nBo,nrec))
    elif recAxis==0:
        zGC = np.zeros((nrec,nBo,2))
        zXY = np.zeros((nrec,nBo,2))
        ztm = np.zeros((nrec,nBo  ))

    jBo = -1
    for jB in range(nBi):
        if pmask[jB] == 1:
            jBo = jBo+1
            zNam[jBo] = pNam[jB]            
            zIDs[jBo] = pIDs[jB]
            for jr in range(nrec):
                if recAxis==2:
                    zGC[jBo,:,jr] =   pGC[jB,:,jr]
                    zXY[jBo,:,jr] =   pXY[jB,:,jr]
                    ztm[jBo,  jr] = ptime[jB  ,jr]
                elif recAxis==0:
                    zGC[jr,jBo,:] =   pGC[jr,jB,:]
                    zXY[jr,jBo,:] =   pXY[jr,jB,:]
                    ztm[jr,jBo  ] = ptime[jr,jB  ]
                    
    return zNam, zIDs, zGC, zXY, ztm




