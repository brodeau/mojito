import numpy as np

from .geomesh import AnglesOfTriangle

class Triangle:

    def __init__( self, nT, xPntID, xCoor, xNbgh ):
        '''
            * nT:     [integer] number of Triangles
            * xPntID: [(nT,3)   array of integers] the 3 point IDs composing the quad, in counter-clockwize
            * xCoor:  [(nT,3,2) array of floats] the 3 [lon,lat] geographic coordinates "      " [degrees]
            * xNbgh:  [(nT,3)   array of integers] the 3 IDs of the 3 neighbor triangles
        '''
        self.nT        = nT
        self.length    = nT
        # ID:    [(nT)     array of integers] ID of the triangle
        self.ID        = np.array([i for i in range(nT)], dtype=int)
        self.pointIDs  = np.array(xPntID, dtype=int)
        self.pointCoor = np.array(xCoor)
        self.neighbors = np.array(xNbgh, dtype=int)

    def angles( self ):
        ''' Returns the (nT,3) array of the 3 angles (from 1st point to 3rd point, so counter-clockwize) '''
        return np.array( [ AnglesOfTriangle( self.pointCoor[i] )  for i in range(self.nT) ] )


class Quadrangle:

    def __init__( self, nQ, xPntID, xCoor ):
        '''
            * nQ:    [integer] number of Quadrangles
            * xPntID:  [(nQ,4)   array of integers] the 4 point IDs composing the quad, in counter-clockwize
            * xCoor: [(nQ,4,2) array of floats] the 4 [lon,lat] geographic coordinates "      " [degrees]
        '''
        self.nQ        = nQ
        self.length    = nQ
        # ID:    [(nT)     array of integers] ID of the triangle
        self.ID        = np.array([i for i in range(nQ)], dtype=int)
        self.pointIDs  = np.array(xPntID, dtype=int)
        self.pointCoor = np.array(xCoor)

    def angles( self ):
        ''' Returns the (nT,4) array of the 4 angles (from 1st point to 4th point, so counter-clockwize) 
        
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
                  o           o
                 A            B

        '''
        xA = np.zeros((self.nQ,4)) - 999.
        for i in range(self.nQ):
            va1 = AnglesOfTriangle( np.array([ self.pointCoor[i,j]  for j in [0,1,3] ]) ); # angles of triangle [ABD]
            va2 = AnglesOfTriangle( np.array([ self.pointCoor[i,j]  for j in [2,3,1] ]) ); # angles of triangle [CDB]
            xA[i,:] = [ va1[0], va1[1]+va2[2], va2[0], va1[2]+va2[1] ]
            #
        return xA

