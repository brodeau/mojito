import numpy as np

from .geomesh import LengthsOfTriangle, AnglesOfTriangle, AreaOfTriangle
from sys import exit


class Triangle:

    def __init__( self,  xPcoor, xTPntIdx, xTnbgh, vPIDs, vPnames ):
        '''
               `nP` points => `nT` Triangles!

            * xPcoor:   [(nP,2] array of floats]  the coordinates of the nP points that define the Triangles
            * xTPntIdx: [(nT,3] array of integers] the 3 point indices composing the triangle, in counter-clockwize
            * xTnbgh:   [(nT,3] array of integers] the 3 indices of the 3 neighbor triangles
            * vPIDs:    [(nP)  vector of integers] an integer to identify each point as in xPcoor
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
        self.PointNames   = np.array( vPnames, dtype='U32' ) ; # point names => shape = (nP)
        self.PointIdx     = zvPntIdx                         ; # indices of all the points making the triangles => shape = (nP)
        self.PointXY      = np.array(xPcoor)                             ; # Coordinates of all the points making the triangles => shape = (nP,2)
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

    def __init__( self,  xPcoor, xQPntIdx, vPIDs, vQnames, date='unknown' ):
        '''
               => `nQ` Quadrangles!

            * xPcoor:   [(nP,2] array of floats] the coordinates of the nP points that define the Quads
            * xQPntIdx: [(nQ,4] array of integers] the 4 point indices composing the quad, in counter-clockwize
            * vPIDs:    [(nP)  vector of integers] an integer to identify each point as in xPcoor
            * vQnames:  [(nQ)  vector of strings]  a string to identify each quadrangle
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

        # Integers:
        self.nP        = nP ; # number of points making the quadrangles
        self.nQ        = nQ ; # number of quadrangles (nQ < nP)
        self.length    = nQ ; #    "        "
        self.date      = date ;  # as a string!

        # NumPy Arrays:
        self.PointIDs     = np.array( vPIDs,   dtype=int )   ; # point IDs   => shape = (nP)        
        self.PointIdx     = zvPntIdx                          ; # indices of all the points making the quadrangles => shape = (nP)
        self.PointXY      = np.array(xPcoor)                             ; # Coordinates of all the points making the quadrangles => shape = (nP,2)
        self.QuaIDs       = np.array([i for i in range(nQ)], dtype=int); # IDs of the quadrangles => shape = (nQ)
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
            zA[i] = np.sum( np.array([  zX[i,      k]*zY[i,(k+1)%4]
                                      - zX[i,(k+1)%4]*zY[i,      k] for k in range(4) ]) )
        del zX,zY
        return 0.5*zA

    def MeshVrtcPntIDs( self ):
        ''' Returns the same as `MeshVrtcPntIdx` but with Point IDs rather than Point indices !'''
        zVpIDs = np.zeros((self.nQ,4), dtype=int)
        for i in range(self.nQ):
            v4pidx = self.MeshVrtcPntIdx[i,:]
            v4pids = self.PointIDs[v4pidx]
            zVpIDs[i,:] = v4pids
        return zVpIDs

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
                             PointIDs=Poly.PointIDs, QuadNames=Poly.QuadNames )
        print('\n *** Quadrangle mesh saved into "'+cfile+'" !')

    if ctype=='T':
        np.savez_compressed( cfile, date=cdate, PointXY=Poly.PointXY, MeshVrtcPntIdx=Poly.MeshVrtcPntIdx,
                             NeighborIDs=Poly.NeighborIDs, PointIDs=Poly.PointIDs, PointNames=Poly.PointNames )
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

    (nPoly,nVrtc) = np.shape(MeshVrtcPntIdx)

    if ctype=='T':
        if nVrtc!=3:
            print('ERROR: [polygons.LoadPolygon()] => wrong number of vertices for a triangle:',nVrtc); exit(0)
        PointNames  = data['PointNames']
        NeighborIDs = data['NeighborIDs'] ; # shape: (nPoly,3)
        POLY = Triangle( PointXY, MeshVrtcPntIdx, NeighborIDs, PointIDs, PointNames )

    if ctype=='Q':
        if nVrtc!=4:
            print('ERROR: [polygons.LoadPolygon()] => wrong number of vertices for a quadrangle:',nVrtc); exit(0)
        QuadNames   = data['QuadNames']
        POLY = Quadrangle( PointXY, MeshVrtcPntIdx, PointIDs, QuadNames, date=cdate )

    return POLY
