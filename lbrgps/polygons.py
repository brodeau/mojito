import numpy as np

from .geomesh import LengthsOfTriangle, AnglesOfTriangle, AreaOfTriangle



class Triangle:

    def __init__( self, nT, xPntID, xCoor, xNbgh ):
        '''
            * nT:     [integer] number of Triangles
            * xPntID: [(nT,3)   array of integers] the 3 point IDs composing the triangle, in counter-clockwize
            * xCoor:  [(nT,3,2) array of floats] the 3 [lon,lat] geographic coordinates "      " [degrees]
            * xNbgh:  [(nT,3)   array of integers] the 3 IDs of the 3 neighbor triangles

                         C
                         o
                        /\
                       /  \
                      /    \
                     /      \
                    /        \
                   /__________\
                  o            o
                 A             B

        '''
        vpIDs, vunqixd = np.unique( xPntID.flatten(), return_index=True )
        vpx   = xCoor[:,:,0].flatten()
        vpy   = xCoor[:,:,1].flatten()
        XY = np.array( [ vpx[vunqixd], vpy[vunqixd] ] ).T
        del vpx, vpy, vunqixd

        # Integers:
        self.nP        = len(vpIDs) ; # number of points making the triangles
        self.nT        = nT ; # number of triangles (nT > nP)
        self.length    = nT ; #    "        "
        #
        # NumPy Arrays:
        self.PointIDs    = vpIDs                          ; # IDs of all the points making the triangles => shape = (nP)
        self.PointXY     = XY                             ; # Coordinates of all the points making the triangles => shape = (nP,2)
        #
        self.TriIDs      = np.array([i for i in range(nT)], dtype=int); # IDs of the triangles => shape = (nT)
        self.MeshPointIDs = np.array(xPntID, dtype=int)                ; # 3 point IDs composing the triangles, CCW => shape = (nT,3)
        self.TriPointXY  = np.array(xCoor)              ; # Coordinates of the 3 points composing the triangle => shape = (nT,3,2)
        self.neighbors   = np.array(xNbgh, dtype=int)
        #
        del vpIDs, XY

        
    def lengths( self ):
        ''' Returns the shape(nT,3) array of the length of the 3 segments defining the triangle (counter-clockwize from 1st point) '''
        return np.array( [ LengthsOfTriangle( self.TriPointXY[i] )  for i in range(self.nT) ] )

    def angles( self ):
        ''' Returns the shape(nT,3) array of the 3 angles (counter-clockwize from 1st point) '''
        return np.array( [ AnglesOfTriangle( self.TriPointXY[i] )  for i in range(self.nT) ] )

    def area( self ):
        ''' Returns the shape(nT) array of the area of the triangles '''
        return np.array( [ AreaOfTriangle( self.TriPointXY[i] )  for i in range(self.nT) ] )
    



class Quadrangle:

    def __init__( self, nQ, xPntID, xCoor ):
        '''
            * nQ:    [integer] number of Quadrangles
            * xPntID:  [(nQ,4)   array of integers] the 4 point IDs composing the quad, in counter-clockwize
            * xCoor: [(nQ,4,2) array of floats] the 4 [lon,lat] geographic coordinates "      " [degrees]


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
        vpIDs, vunqixd = np.unique( xPntID.flatten(), return_index=True )
        nP = len(vpIDs) ; # number of points that defines the nQ quadrangles
        
        vpx   = xCoor[:,:,0].flatten()
        vpy   = xCoor[:,:,1].flatten()
        XY = np.array( [ vpx[vunqixd], vpy[vunqixd] ] ).T
        del vpx, vpy, vunqixd

        #print('LOLO: xPntID =', xPntID)
        #print('LOLO:  vpIDs =', vpIDs, nP)

        zPntIdx = xPntID.copy(); # La dedans, je veux foutre les nouvelles IDs (allant de 0 a nP==len(vpIDs), qui est plus petit que le nombres de points initial!!!)
        zPntIdx[:,:] = -1
        for jP in range(nP):
            ii = vpIDs[jP]
            idx = np.where(xPntID==ii)
            zPntIdx[idx] = jP
        #print('LOLO: zPntIdx =', zPntIdx) ; exit(0)
        
        # Integers:
        self.nP        = nP ; # number of points making the quadrangles
        self.nQ        = nQ ; # number of quadrangles (nT > nP)        
        self.length    = nQ ; #    "        "

        # NumPy Arrays:
        self.PointIDs    = vpIDs                          ; # IDs of all the points making the quadrangles => shape = (nP)
        self.PointXY     = XY                             ; # Coordinates of all the points making the quadrangles => shape = (nP,2)
        #
        self.QuaIDs      = np.array([i for i in range(nQ)], dtype=int); # IDs of the quadrangles => shape = (nQ)
        self.MeshPointIDs = np.array(xPntID, dtype=int)                ; # 4 point IDs composing the quadrangles, CCW => shape = (nQ,4)
        self.MeshPointTriIDs = zPntIdx                                    ; # 4 point IDs composing the quadrangles BUT with vpIDs !!!!
        self.QuaPointXY  = np.array(xCoor)              ; # Coordinates of the 4 points composing the quadrangle => shape = (nQ,4,2)

        del vpIDs, XY, zPntIdx

    def lengths( self ):
        ''' Returns the shape(nQ,4) array of the length of the 4 segments defining the triangle (counter-clockwize from 1st point) '''
        zL = np.zeros((self.nQ,4)) - 999.
        for i in range(self.nQ):
            vL1 = LengthsOfTriangle( np.array([ self.QuaPointXY[i,j]  for j in [0,1,3] ]) ); # length of sides of triangle [ABD]
            vL2 = LengthsOfTriangle( np.array([ self.QuaPointXY[i,j]  for j in [2,3,1] ]) ); # length of sides of triangle [CDB]
            zL[i,:] = [ vL1[0], vL2[2], vL2[0], vL1[2] ] ; # lengths of [AB], [BC], [CD], [DA]
            #
        return zL
        
    def angles( self ):
        ''' Returns the shape(nQ,4) array of the 4 angles (from 1st point to 4th point, so counter-clockwize)         
        '''
        za = np.zeros((self.nQ,4)) - 999.
        for i in range(self.nQ):
            va1 = AnglesOfTriangle( np.array([ self.QuaPointXY[i,j]  for j in [0,1,3] ]) ); # angles of triangle [ABD]
            va2 = AnglesOfTriangle( np.array([ self.QuaPointXY[i,j]  for j in [2,3,1] ]) ); # angles of triangle [CDB]
            za[i,:] = [ va1[0], va1[1]+va2[2], va2[0], va1[2]+va2[1] ] ; # angles at A, B, C, D
            #
        return za

    def area( self ):
        ''' Returns the shape(nQ) array of the area of the quadrangles '''
        zA = np.zeros(self.nQ) - 999.
        for i in range(self.nQ):
            ra1 = AreaOfTriangle( np.array([ self.QuaPointXY[i,j]  for j in [0,1,3] ]) ); # area of triangle [ABD]
            ra2 = AreaOfTriangle( np.array([ self.QuaPointXY[i,j]  for j in [2,3,1] ]) ); # area of triangle [CDB]
            zA[i] = ra1 + ra2
            #
        return zA

        
    #return np.array( [ AreaOfTriangle( self.TriPointXY[i] )  for i in range(self.nT) ] )






def SavePolygon( cfile, Poly, ctype='Q' ):
    '''
        cfile: file to save into
        Poly:  polygon object to save (`polygon.Triangle` or `polygon.Quadrangle`)
        ctype: type of polygon object: 'Q' => Quadrangle, 'T' => Triangle

    '''
    if not ctype in ['Q','T']:
        print('ERROR: [polygons.SavePolygon()] => wrong polygon type'); exit(0)

    if ctype=='Q':
        np.savez_compressed( cfile, PointXY=Poly.PointXY, Mesh=Poly.MeshPointIDs, MeshTriPntIDs=Poly.MeshPointTriIDs,
                             Lengths=Poly.lengths(), Angles=Poly.angles(), Areas=Poly.area() )
        print('\n *** Quadrangle mesh saved into "'+cfile+'" !')

    if ctype=='T':
        np.savez_compressed( cfile, PointXY=Poly.PointXY, Mesh=Poly.MeshPointIDs,
                             Lengths=Poly.lengths(), Angles=Poly.angles(), Areas=Poly.area() ) ; #, names=vnam )
        print('\n *** Triangle mesh saved into "'+cfile+'" !')


# Loading:
# Qxy    => PointXY
# Qmesh  => Mesh
# Qmeix  => MeshTriPntIDs
#
# Txy    => PointXY
# Tmesh  => Mesh

# Class:
# QuaPointIDs => MeshPointIDs
# QuaPointIdx => MeshPointTriIDs
#
# TriPointIDs => MeshPointIDs
