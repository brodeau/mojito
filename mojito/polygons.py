import numpy as np

from .geomesh import LengthsOfTriangle, AnglesOfTriangle, AreaOfTriangle, AreaOfQuadrangle
from sys import exit


class Triangle:

    def __init__( self,  xPcoor, xTPntIdx, xTnbgh, vPIDs, vPtime, vPnames ):
        '''
               `nP` points => `nT` Triangles!

            * xPcoor:   [(nP,2] array of floats]  the coordinates of the nP points that define the Triangles
            * xTPntIdx: [(nT,3] array of integers] the 3 point indices composing the triangle, in counter-clockwize
            * xTnbgh:   [(nT,3] array of integers] the 3 indices of the 3 neighbor triangles
            * vPIDs:    [(nP)  vector of integers] an integer to identify each point as in xPcoor
            * vPtime:   [(nP)  vector of floats] epoch/unix time corresponding to buoy's position
            * vPnames:  [(nP)  vector of strings]  a string to identify each point as in xPcoor

                         C
                         o
                        /\
                       /  \
                      /    \
                     /      \
                    /        \
                   /__________\
                  o            o
                  A            B

        '''
        cEM = 'ERROR: [polygons.Triangle] =>'
        #
        (nP,nd2) = np.shape(xPcoor)
        if nd2!=2:
            print(cEM+' problem in the shape of coordinate array!'); exit(0)
        (nT,nd3) = np.shape(xTPntIdx)
        if nd3!=3:
            print(cEM+' we expected 2nd axis to have length 3, not:',nd3,' in the Triangle array!'); exit(0)
        if np.shape(xTnbgh) != (nT,3):
            print(cEM+' problem in the shape of neighbor array!'); exit(0)
        if (nP,) != np.shape(vPnames):
            print(cEM+' problem in the length of point names!'); exit(0)
            
        zvPntIdx = np.unique( xTPntIdx.flatten() ) ; # Populate all the points (IDs) involved
        if len(zvPntIdx) != nP:
            print(cEM+' problem with the number of points at play deduced from `xTPntIdx`', len(zvPntIdx), nP); exit(0)

        zTcoor = np.array([ [ xPcoor[i,:] for i in xTPntIdx[jT,:] ] for jT in range(nT) ]) ; # for each Triangle the 3 coordinates of 3 points [nQ,3,2]

        # Integers:
        self.nP        = nP ; # number of points making the triangles
        self.nT        = nT ; # number of triangles (nT < nP)
        self.length    = nT ; #    "        "
        self.date      = 'unknown' ;  # we are not going to work with triangles...
                
        # NumPy Arrays:
        self.PointIDs     = np.array( vPIDs,   dtype=int )   ; # point IDs   => shape = (nP)
        self.PointTime    = np.array( vPtime )               ; # point time  => shape = (nP)
        self.PointNames   = np.array( vPnames, dtype='U32' ) ; # point names => shape = (nP)
        self.PointIdx     = zvPntIdx                         ; # indices of all the points making the triangles => shape = (nP)
        self.PointXY      = np.array(xPcoor)                 ; # Coordinates of all the points making the triangles => shape = (nP,2)
        self.TriIDs       = np.array([i for i in range(nT)], dtype=int); # IDs of the triangles => shape = (nT)
        self.MeshVrtcPntIdx = np.array(xTPntIdx, dtype=int)                ; # 3 point indices composing the triangles, CCW => shape = (nT,3)
        self.MeshPointXY  = np.array(zTcoor)              ; # Coordinates of the 3 points composing the triangle => shape = (nT,3,2)
        self.NeighborIDs  = np.array(xTnbgh, dtype=int)   ; # IDs of neighbor triangles...
        #
        del zvPntIdx, zTcoor

    def lengths( self ):
        ''' Returns [nT,3] array of the length of the 3 segments defining the triangle (counter-clockwize from 1st point) '''
        return np.array( [ LengthsOfTriangle( self.MeshPointXY[i] )  for i in range(self.nT) ] )

    def angles( self ):
        ''' Returns [nT,3] array of the 3 angles (counter-clockwize from 1st point) '''
        return np.array( [ AnglesOfTriangle( self.MeshPointXY[i] )  for i in range(self.nT) ] )

    def area( self ):
        ''' Returns [nT] array of the area of the triangles '''
        return np.array( [ AreaOfTriangle( self.MeshPointXY[i] )  for i in range(self.nT) ] )




class Quadrangle:

    def __init__( self,  xPcoor, xQPntIdx, vPIDs, vPtime, vQnames, vQIDs=[], date='unknown' ):
        '''
               => `nQ` Quadrangles!

            * xPcoor:   [(nP,2] array of floats] the coordinates of the nP points that define the Quads
            * xQPntIdx: [(nQ,4] array of integers] the 4 point indices composing the quad, in counter-clockwize
            * vPIDs:    [(nP)  vector of integers] an integer to identify each point as in xPcoor
            * vPtime:    [(nP)  vector of floats] epoch/unix time corresponding to buoy's position
            * vQnames:  [(nQ)  vector of strings]  a string to identify each quadrangle
            * vQIDs:    [(nQ)  vector of integers] force IDs of quadrangle to this instead of using 0 to nQ !!!
            * date:   OPTIONAL string of the form `YYYY-MM-DD_hh:mm:ss`

            IMPORTANT: we expect the 1st and 3rd (indices 0 & 2) elements of the 4 points to be the
                       the two points of the quadrangle facing the diagonal that was the common segment
                       between the 2 triangles that formed this quarangle

                         D           C
                         o___________o
                        /\          /
                       /  \        /
                      /    \      /
                     /      \    /     => here it would be A and C, the common segment being B-D
                    /        \  /
                   /__________\/
                  o            o
                  A            B

        '''
        cEM = 'ERROR: [polygons.Quadrangle] =>'
        #        
        (nP,nd2) = np.shape(xPcoor)
        if nd2!=2:
            print(cEM+' problem in the shape of coordinate array!'); exit(0)
        (nQ,nd4) = np.shape(xQPntIdx)
        if nd4!=4:
            print(cEM+' we expected 2nd axis to have length 4, not:',nd4,' in the Quad array!'); exit(0)
        if (nQ,) != np.shape(vQnames):
            print(cEM+' problem in the length of point names!', nP, np.shape(vQnames)); exit(0)

        zvPntIdx = np.unique( xQPntIdx.flatten() ) ; # Populate all the points (IDs) involved
        if len(zvPntIdx) != nP:
            print(cEM+' problem with the number of points at play deduced from `xQPntIdx`: ',len(zvPntIdx),nP); exit(0)

        zTcoor = np.array([ [ xPcoor[i,:] for i in xQPntIdx[jQ,:] ] for jQ in range(nQ) ]) ; # for each Quadrangle the 4 coordinates of 4 points [nQ,4,2]

        l_force_Q_IDs = (len(vQIDs) == nQ)

        # Integers:
        self.nP        = nP ; # number of points making the quadrangles
        self.nQ        = nQ ; # number of quadrangles (nQ < nP)
        self.length    = nQ ; #    "        "
        self.date      = date ;  # as a string!

        # NumPy Arrays:
        self.PointIDs     = np.array( vPIDs,   dtype=int )   ; # point IDs   => shape = (nP)
        self.PointTime    = np.array( vPtime )               ; # point time  => shape = (nP)        
        self.PointIdx     = zvPntIdx                         ; # indices of all the points making the quadrangles => shape = (nP)
        self.PointXY      = np.array(xPcoor)                             ; # Coordinates of all the points making the quadrangles => shape = (nP,2)
        if l_force_Q_IDs:
            self.QuadIDs       = np.array(vQIDs, dtype=int); # IDs of the quadrangles => shape = (nQ)
        else:
            self.QuadIDs       = np.array([i for i in range(nQ)], dtype=int); # IDs of the quadrangles => shape = (nQ)
        self.MeshVrtcPntIdx = np.array(xQPntIdx, dtype=int)   ; # 4 point indices composing the quadrangles (indices in the frame of "only Quadrangles")
        self.MeshPointXY  = np.array(zTcoor)              ; # Coordinates of the 4 points composing the quadrangle => shape = (nQ,4,2)
        self.QuadNames   = np.array( vQnames, dtype='U32' ) ; # point names => shape = (nP)
        #
        del zvPntIdx, zTcoor


    def lengths( self ):
        ''' Returns [nQ,4] array of the length of the 4 segments defining the quadrangle (counter-clockwize from 1st point) '''
        zL = np.zeros((self.nQ,4)) - 999.
        for i in range(self.nQ):
            vL1 = LengthsOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [0,1,3] ]) ); # length of sides of triangle [ABD]
            vL2 = LengthsOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [2,3,1] ]) ); # length of sides of triangle [CDB]
            zL[i,:] = [ vL1[0], vL2[2], vL2[0], vL1[2] ] ; # lengths of [AB], [BC], [CD], [DA]
            #
        return zL

    def angles( self ):
        ''' Returns [nQ,4] array of the 4 angles (from 1st point to 4th point, so counter-clockwize)
        '''
        za = np.zeros((self.nQ,4)) - 999.
        for i in range(self.nQ):
            va1 = AnglesOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [0,1,3] ]) ); # angles of triangle [ABD]
            va2 = AnglesOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [2,3,1] ]) ); # angles of triangle [CDB]
            za[i,:] = [ va1[0], va1[1]+va2[2], va2[0], va1[2]+va2[1] ] ; # angles at A, B, C, D
            #
        return za

    def area( self ):
        ''' Returns [nQ] array of the area of the quadrangles, see for example Bouchat et al. 2022 (Eq.5)'''
        zX , zY = self.MeshPointXY[:,:,0].copy() , self.MeshPointXY[:,:,1].copy()
        zA      = np.zeros(self.nQ) - 9999.
        #
        for i in range(self.nQ):
            zA[i] = AreaOfQuadrangle( zX[i,:] , zY[i,:])
            #zA[i] = np.sum( np.array([  zX[i,      k]*zY[i,(k+1)%4]
            #                          - zX[i,(k+1)%4]*zY[i,      k] for k in range(4) ]) )
        del zX,zY
        #return 0.5*zA
        return zA

    def MeshVrtcPntIDs( self ):
        ''' Returns the same as `MeshVrtcPntIdx` but with Point IDs rather than Point indices !'''
        zVpIDs = np.zeros((self.nQ,4), dtype=int)
        for i in range(self.nQ):
            v4pidx = self.MeshVrtcPntIdx[i,:]
            v4pids = self.PointIDs[v4pidx]
            zVpIDs[i,:] = v4pids
        return zVpIDs

    def MeshVrtcPntTime( self ):
        ''' Returns the same as `MeshVrtcPntIdx` but with Point Time rather than Point indices !'''
        zVpTime = np.zeros((self.nQ,4))
        for i in range(self.nQ):
            v4pidx = self.MeshVrtcPntIdx[i,:]
            v4pids = self.PointTime[v4pidx]
            zVpTime[i,:] = v4pids
        return zVpTime

###################################################################################################################




def SaveClassPolygon( cfile, Poly, ctype='Q', force_date='unknown' ):
    '''
        Save all arrays necessary to rebuild the Polygon object later on.

       * cfile: file to save into
       * Poly:  polygon object to save (`polygon.Triangle` or `polygon.Quadrangle`)
       * ctype: type of polygon object: 'Q' => Quadrangle, 'T' => Triangle
       * force_date:   OPTIONAL string of the form `YYYY-MM-DD_hh:mm:ss`

    '''
    if not ctype in ['Q','T']:
        print('ERROR: [polygons.SavePolygon()] => wrong polygon type'); exit(0)

    cdate = Poly.date
    if cdate != 'unknown' and force_date != 'unknown':
        print('WARNING: [polygons.SavePolygon()] => overriding class original date: '+cdate+' with: '+force_date); exit(0)
        cdate = force_date
    if ctype=='Q':
        np.savez_compressed( cfile, date=cdate, PointXY=Poly.PointXY, MeshVrtcPntIdx=Poly.MeshVrtcPntIdx,
                             PointIDs=Poly.PointIDs, PointTime=Poly.PointTime, QuadNames=Poly.QuadNames, QuadIDs=Poly.QuadIDs )
        print('\n *** Quadrangle mesh saved into "'+cfile+'" !')

    if ctype=='T':
        np.savez_compressed( cfile, date=cdate, PointXY=Poly.PointXY, MeshVrtcPntIdx=Poly.MeshVrtcPntIdx, NeighborIDs=Poly.NeighborIDs,
                             PointIDs=Poly.PointIDs, PointTime=Poly.PointTime, PointNames=Poly.PointNames )
        print('\n *** Triangle mesh saved into "'+cfile+'" !')




def LoadClassPolygon( cfile, ctype='Q' ):
    '''
        Recreate the polygon object out of what's read into the `npz` (saved via `SaveClassPolygon()`)

        cfile: file to save into
        Poly:  polygon object to save (`polygon.Triangle` or `polygon.Quadrangle`)
        ctype: type of polygon object: 'Q' => Quadrangle, 'T' => Triangle

    '''
    if not ctype in ['Q','T']:
        print('ERROR: [polygons.LoadPolygon()] => wrong polygon type'); exit(0)


    print('\n\n *** Loading the "'+ctype+'"class from "'+cfile+'" ...')

    data = np.load(cfile) ; #, allow_pickle=True)

    cdate          = data['date']
    PointXY        = data['PointXY']
    MeshVrtcPntIdx = data['MeshVrtcPntIdx']        ; # the `nVrtc` point indices for each polygon => shape: (nPoly,nVrtc)
    PointIDs       = data['PointIDs']
    PointTime      = data['PointTime']

    (nPoly,nVrtc) = np.shape(MeshVrtcPntIdx)

    if ctype=='T':
        if nVrtc!=3:
            print('ERROR: [polygons.LoadPolygon()] => wrong number of vertices for a triangle:',nVrtc); exit(0)
        PointNames  = data['PointNames']
        NeighborIDs = data['NeighborIDs'] ; # shape: (nPoly,3)
        POLY = Triangle( PointXY, MeshVrtcPntIdx, NeighborIDs, PointIDs, PointTime, PointNames )

    if ctype=='Q':
        if nVrtc!=4:
            print('ERROR: [polygons.LoadPolygon()] => wrong number of vertices for a quadrangle:',nVrtc); exit(0)
        QuadNames = data['QuadNames']
        QuadIDs   = data['QuadIDs']
        POLY = Quadrangle( PointXY, MeshVrtcPntIdx, PointIDs, PointTime, QuadNames, vQIDs=QuadIDs, date=cdate )

    return POLY



def RecycleQuads( pCoor, pTime, pIDs, pQDS,  iverbose=0 ):
    '''
          pCoor:  new coordinates of the cloud of points (nP,2)
          pTime:  new time for these points              (nP)
          pIDs :  IDs of points as in `pCoor`            (nP)
          pQDS :  Quad class based on the previous cloud of points

        => Returns: all arrays necessary to define the new Quad class based on the new location of the points
    '''
    cEM = 'ERROR: [polygons.RecycleQuads] =>'
    
    nP0     = pQDS.nP
    nQ0     = pQDS.nQ

    zQVids_0 = pQDS.MeshVrtcPntIdx.copy() ; # (nQ,4)
    zPids_0  = pQDS.PointIDs.copy()       ; # (nP) All point IDs involved in initial polpulation of Quads
    zQnms_0  = pQDS.QuadNames.copy()      ; # (nQ)
    zQids_0  = pQDS.QuadIDs.copy()        ; # (nQ)   All Quad IDs involved in initial polpulation of Quads

    lst_PIDs_aliv = []
    lst_QIDs_aliv = []

    for jQ in range(nQ0):
        l_Survived = True
        idQuad = zQids_0[jQ]
        print('\n --- Quad #'+str(jQ)+' with ID '+str(idQuad)+':')
    
        # Indices & IDs of the 4 points involved in this Quad:
        v4Pind = zQVids_0[jQ,:]
        v4Pids = np.array( [ zPids_0[i] for i in v4Pind ], dtype=int )

        # If any of the 4 points involved in this Quad does not exist in current cloud point:
        for jid in v4Pids:
            if not jid in pIDs:
                print('     => point with ID '+str(jid)+' does not exist anymore in this record')
                print('     ===> quad with ID '+str(idQuad)+' is forgotten...')
                l_Survived = False
                break
    
        if l_Survived:
            # Current Quad still exists in this record...
            if iverbose>0: print('     => still exists at current record :)!')
            lst_QIDs_aliv.append(idQuad)
            lst_PIDs_aliv.append(v4Pids)

    zQids_aliv = np.array(lst_QIDs_aliv)
    zPIDs_aliv = np.unique(np.array(lst_PIDs_aliv)) ; # we probably break the initial order...
    nQ = len(zQids_aliv)
    nP = len(zPIDs_aliv)

    if iverbose>1:
        print('  => IDs of  quads still involved:', zQids_aliv)                
        print('  => IDs of points still involved:', zPIDs_aliv)                
    del v4Pind, v4Pids, lst_QIDs_aliv, lst_PIDs_aliv
    
    print('\n *** At present record, we have '+str(nQ)+' Quads / '+str(nQ0)+' that survived!')
    print('        and  '+str(nP)+' points  / '+str(nP0)+' involved...')
    
    # * xQpnts remains the same! That's the whole point!!!
    # * xQcoor should be updated with the new point coordinates at this record:

    if nQ < nQ0:
        print('\n *** nQ < nQ0 => need to shrink arrays! :(')
        zPcoor = np.zeros((nP,2))              ; # coordinates of the nP points
        zPtime = np.zeros( nP   )              ; # time of the nP points
        zQVids = np.zeros((nQ,4), dtype=int)   ; # indices (to be used in xCoor) of the 4 points composing the Quad
        zPids  = np.zeros(nP)                  ; # IDs of points
        zQnms  = np.zeros(nQ,     dtype='U32') ; # Quad names
        zQids  = np.zeros(nQ,     dtype=int)   ; # Quad IDs

        # Loop along points:
        lst_pntID_used = []
        jP=-1
        for jPp in range(nP0):
            idPoint = zPids_0[jPp]
            if idPoint in zPIDs_aliv and not idPoint in lst_pntID_used:
                jP = jP+1
                lst_pntID_used.append(idPoint)
                zPids[jP] = idPoint

                (ip,) = np.where(pIDs==idPoint)
                if len(ip)!=1: print(cEM+'  Z1'); exit(0)
                zPcoor[jP,:] = pCoor[ip[0],:]
                zPtime[jP]   = pTime[ip[0]]
        del lst_pntID_used
        if jP != nP-1: print(cEM+'  Z2'); exit(0)

        # Loop along quads:
        jQ=-1
        for jQp in range(nQ0):
            #
            l_Survived = True
            idQuad = zQids_0[jQp]
            # Indices & IDs of the 4 points involved in this Quad:
            v4Pind = zQVids_0[jQp,:]
            v4Pids = np.array( [ zPids_0[i] for i in v4Pind ], dtype=int )
            for jid in v4Pids:
                if not jid in pIDs:
                    l_Survived = False
                    break
            if l_Survived:
                # Current Quad still exists in this record...
                jQ = jQ+1
                if not idQuad in zQids_aliv: print(cEM+'  Y2'); exit(0)                                                        
                zQnms[jQ]    = zQnms_0[jQp]
                zQids[jQ]    = zQids_0[jQp]
                # Tricky for `xQpnts` as we must also update the new point indices in the new frame [0,nP-1]
                for jp in range(4):
                    iPidx = v4Pind[jp]         ; # this index is valid for 0ious go, i.e. in the frame [0,nP0-1]
                    iPids = zPids_0[iPidx]  ; # => IDs of the corresponding point
                    (ip,) = np.where( zPids == iPids ) ; # => index of this point in the frame [0,nP-1] ; #fixme
                    if len(ip)!=1: print(cEM+'  Z5'); exit(0)
                    zQVids[jQ,jp] = ip[0]

                    
    elif nQ == nQ0:
        print('\n *** nQ == nQ0 => NO need to shrink arrays! :D ')
        
        # Only zPcoor must be modified, taking into consideration new coordinates of the point cloud        
        zQVids = zQVids_0
        zPids = zPids_0
        zQnms  = zQnms_0
        zQids  = zQids_0
        #
        zPcoor = np.zeros((nP,2))
        zPtime = np.zeros( nP )
        _,ind2keep,_ = np.intersect1d(pIDs, zPids, return_indices=True); # retain only indices of `pIDs` that exist in `zPids`
        zPcoor[:,:] = pCoor[ind2keep,:]
        zPtime[:]   = pTime[ind2keep]
        
    else:
        print(cEM+' : nQ > nQ0!!!'); exit(0)

    return zPcoor, zPtime, zQVids, zPids, zQnms, zQids

