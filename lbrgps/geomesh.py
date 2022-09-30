import numpy as nmp

class Triangle:

    def __init__( self, nT, xPntID, xCoor, xNbgh ):
        '''
            * nT:     [integer] number of Triangles
            * xPntID: [(nT,3)   array of integers] the 3 point IDs composing the quad, in counter-clockwize
            * xCoor:  [(nT,3,2) array of floats] the 3 [lon,lat] geographic coordinates "      " [degrees]
            * xNbgh:  [(nT,3)   array of integers] the 3 IDs of the 3 neighbor triangles
        '''
        self.length    = nT
        # ID:    [(nT)     array of integers] ID of the triangle
        self.ID        = nmp.array([i for i in range(nT)])
        self.pointIDs  = nmp.array(xPntID, dtype=int)
        self.pointCoor = nmp.array(xCoor)
        self.neighbors = nmp.array(xNbgh, dtype=int)

    #def angles( self ):
        # returns a (nQ,4) arrays of the 4 angles (from 1st point to 4th point, so clockwize)




class Quadrangle:

    def __init__( self, nQ, xPntID, xCoor ):
        '''
            * nQ:    [integer] number of Quadrangles
            * xPntID:  [(nQ,4)   array of integers] the 4 point IDs composing the quad, in counter-clockwize
            * xCoor: [(nQ,4,2) array of floats] the 4 [lon,lat] geographic coordinates "      " [degrees]
        '''
        self.length    = nQ
        # ID:    [(nT)     array of integers] ID of the triangle
        self.ID        = nmp.array([i for i in range(nQ)])
        self.pointIDs  = nmp.array(xPntID, dtype=int)
        self.pointCoor = nmp.array(xCoor)

    #def angles( self ):
        # returns a (nQ,4) arrays of the 4 angles (from 1st point to 4th point, so clockwize)

