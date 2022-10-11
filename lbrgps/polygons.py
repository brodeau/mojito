import numpy as np

from .geomesh import LengthsOfTriangle, AnglesOfTriangle, AreaOfTriangle
from sys import exit


### TO DO:
###
### what we presently refer to "ID" should become "index/indices", the real IDs should be in `cNames`...


class Triangle:

    def __init__( self,  xPcoor, xTpntID, xTnbgh, vPnames ):
        '''
               `nP` points => `nT` Triangles!

            * xPcoor:  [(nP,2) array of floats]  the coordinates of all points that define the Quads
            * xTpntID: [(nT,3) array of integers] the 3 point IDs composing the triangle, in counter-clockwize
            * xTnbgh:  [(nT,3) array of integers] the 3 IDs of the 3 neighbor triangles
            * vPnames: [(nP)  vector of strings]  a string to identify each point

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
        (nP,nd2) = np.shape(xPcoor)
        if nd2!=2:
            print('ERROR: [polygons.Triangle] => problem in the shape of coordinate array!'); exit(0)
        (nT,nd3) = np.shape(xTpntID)
        if nd3!=3:
            print('ERROR: [polygons.Triangle] => we expected 2nd axis to have length 3, not:',nd3,' in the Triangle array!'); exit(0)
        if np.shape(xTnbgh) != (nT,3):
            print('ERROR: [polygons.Triangle] => problem in the shape of neighbor array!'); exit(0)
        if (nP,) != np.shape(vPnames):
            print('ERROR: [polygons.Triangle] => problem in the length of point names!'); exit(0)
            
        zvIDs = np.unique( xTpntID.flatten() )
        if len(zvIDs) != nP:
            print('ERROR: [polygons.Triangle] => problem with the number of points at play deduced from `xTpntID`'); exit(0)

        zTcoor = np.array([ [ xPcoor[i,:] for i in xTpntID[jT,:] ] for jT in range(nT) ]) ; # for each Triangle the 3 coordinates of 3 points [nQ,3,2]

        # Integers:
        self.nP        = nP ; # number of points making the triangles
        self.nT        = nT ; # number of triangles (nT < nP)
        self.length    = nT ; #    "        "
        # NumPy Arrays:
        self.PointNames   = np.array( vPnames, dtype='U32' ) ; # point names => shape = (nP)
        self.PointIDs     = zvIDs                          ; # IDs of all the points making the triangles => shape = (nP)
        self.PointXY      = np.array(xPcoor)                             ; # Coordinates of all the points making the triangles => shape = (nP,2)
        self.TriIDs       = np.array([i for i in range(nT)], dtype=int); # IDs of the triangles => shape = (nT)
        self.MeshPointIDs = np.array(xTpntID, dtype=int)                ; # 3 point IDs composing the triangles, CCW => shape = (nT,3)
        self.MeshPointXY  = np.array(zTcoor)              ; # Coordinates of the 3 points composing the triangle => shape = (nT,3,2)
        self.NeighborIDs  = np.array(xTnbgh, dtype=int)
        #
        del zvIDs, zTcoor

    def lengths( self ):
        ''' Returns the shape(nT,3) array of the length of the 3 segments defining the triangle (counter-clockwize from 1st point) '''
        return np.array( [ LengthsOfTriangle( self.MeshPointXY[i] )  for i in range(self.nT) ] )

    def angles( self ):
        ''' Returns the shape(nT,3) array of the 3 angles (counter-clockwize from 1st point) '''
        return np.array( [ AnglesOfTriangle( self.MeshPointXY[i] )  for i in range(self.nT) ] )

    def area( self ):
        ''' Returns the shape(nT) array of the area of the triangles '''
        return np.array( [ AreaOfTriangle( self.MeshPointXY[i] )  for i in range(self.nT) ] )




class Quadrangle:

    def __init__( self,  xPcoor, xQpntID, vQnames ):
        '''
               => `nQ` Quadrangles!

            * xPcoor:  [(nP,2) array of floats] the coordinates of all points that define the Quads
            * xQpntID: [(nQ,4) array of integers] the 4 point IDs composing the quad, in counter-clockwize
            * vQnames: [(nQ)  vector of strings]  a string to identify each quadrangle

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
        (nP,nd2) = np.shape(xPcoor)
        if nd2!=2:
            print('ERROR: [polygons.Quadrangle] => problem in the shape of coordinate array!'); exit(0)
        (nQ,nd4) = np.shape(xQpntID)
        if nd4!=4:
            print('ERROR: [polygons.Quadrangle] => we expected 2nd axis to have length 4, not:',nd4,' in the Quad array!'); exit(0)
        if (nQ,) != np.shape(vQnames):
            print('ERROR: [polygons.Quadrangle] => problem in the length of point names!', nP, np.shape(vQnames)); exit(0)

        zvIDs = np.unique( xQpntID.flatten() )
        if len(zvIDs) != nP:
            print('ERROR: [polygons.Quadrangle] => problem with the number of points at play deduced from `xQpntID`'); exit(0)

        zTcoor = np.array([ [ xPcoor[i,:] for i in xQpntID[jQ,:] ] for jQ in range(nQ) ]) ; # for each Quadrangle the 4 coordinates of 4 points [nQ,4,2]

        # Integers:
        self.nP        = nP ; # number of points making the quadrangles
        self.nQ        = nQ ; # number of quadrangles (nQ < nP)
        self.length    = nQ ; #    "        "

        # NumPy Arrays:
        self.PointNames   = np.array( vQnames, dtype='U32' ) ; # point names => shape = (nP)
        self.PointIDs     = zvIDs                          ; # IDs of all the points making the quadrangles => shape = (nP)
        self.PointXY      = np.array(xPcoor)                             ; # Coordinates of all the points making the quadrangles => shape = (nP,2)
        self.QuaIDs       = np.array([i for i in range(nQ)], dtype=int); # IDs of the quadrangles => shape = (nQ)
        self.MeshPointIDs = np.array(xQpntID, dtype=int)   ; # 4 point IDs composing the quadrangles (IDs in the frame of "only Quadrangles")
        self.MeshPointXY  = np.array(zTcoor)              ; # Coordinates of the 4 points composing the quadrangle => shape = (nQ,4,2)
        #
        del zvIDs, zTcoor


    def lengths( self ):
        ''' Returns the shape(nQ,4) array of the length of the 4 segments defining the quadrangle (counter-clockwize from 1st point) '''
        zL = np.zeros((self.nQ,4)) - 999.
        for i in range(self.nQ):
            vL1 = LengthsOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [0,1,3] ]) ); # length of sides of triangle [ABD]
            vL2 = LengthsOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [2,3,1] ]) ); # length of sides of triangle [CDB]
            zL[i,:] = [ vL1[0], vL2[2], vL2[0], vL1[2] ] ; # lengths of [AB], [BC], [CD], [DA]
            #
        return zL

    def angles( self ):
        ''' Returns the shape(nQ,4) array of the 4 angles (from 1st point to 4th point, so counter-clockwize)
        '''
        za = np.zeros((self.nQ,4)) - 999.
        for i in range(self.nQ):
            va1 = AnglesOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [0,1,3] ]) ); # angles of triangle [ABD]
            va2 = AnglesOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [2,3,1] ]) ); # angles of triangle [CDB]
            za[i,:] = [ va1[0], va1[1]+va2[2], va2[0], va1[2]+va2[1] ] ; # angles at A, B, C, D
            #
        return za

    def area( self ):
        ''' Returns the shape(nQ) array of the area of the quadrangles '''
        zA = np.zeros(self.nQ) - 999.
        for i in range(self.nQ):
            ra1 = AreaOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [0,1,3] ]) ); # area of triangle [ABD]
            ra2 = AreaOfTriangle( np.array([ self.MeshPointXY[i,j]  for j in [2,3,1] ]) ); # area of triangle [CDB]
            zA[i] = ra1 + ra2
            #
        return zA



###################################################################################################################




def SaveClassPolygon( cfile, Poly, ctype='Q' ):
    '''
        Save all arrays necessary to rebuild the Polygon object later on.

        cfile: file to save into
        Poly:  polygon object to save (`polygon.Triangle` or `polygon.Quadrangle`)
        ctype: type of polygon object: 'Q' => Quadrangle, 'T' => Triangle

    '''
    if not ctype in ['Q','T']:
        print('ERROR: [polygons.SavePolygon()] => wrong polygon type'); exit(0)

    if ctype=='Q':
        np.savez_compressed( cfile, PointXY=Poly.PointXY, MeshPointIDs=Poly.MeshPointIDs,
                             PointNames=Poly.PointNames )
        print('\n *** Quadrangle mesh saved into "'+cfile+'" !')

    if ctype=='T':
        np.savez_compressed( cfile, PointXY=Poly.PointXY, MeshPointIDs=Poly.MeshPointIDs,
                             NeighborIDs=Poly.NeighborIDs, PointNames=Poly.PointNames )
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

    PointXY      = data['PointXY']
    MeshPointIDs = data['MeshPointIDs']        ; # the `nVrtc` point IDs for each polygon => shape: (nPoly,nVrtc)
    PointNames   = data['PointNames']

    (nPoly,nVrtc) = np.shape(MeshPointIDs)

    if ctype=='T':
        if nVrtc!=3:
            print('ERROR: [polygons.LoadPolygon()] => wrong number of vertices for a triangle:',nVrtc); exit(0)
        NeighborIDs = data['NeighborIDs'] ; # shape: (nPoly,3)
        POLY = Triangle( PointXY, MeshPointIDs, NeighborIDs, PointNames )

    if ctype=='Q':
        if nVrtc!=4:
            print('ERROR: [polygons.LoadPolygon()] => wrong number of vertices for a quadrangle:',nVrtc); exit(0)
        POLY = Quadrangle( PointXY, MeshPointIDs, PointNames )

    return POLY
