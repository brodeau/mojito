import numpy as np

from .geomesh import AnglesOfTriangle



class Triangle:

    def __init__( self, nT, xPntID, xCoor, xNbgh ):
        '''
            * nT:     [integer] number of Triangles
            * xPntID: [(nT,3)   array of integers] the 3 point IDs composing the triangle, in counter-clockwize
            * xCoor:  [(nT,3,2) array of floats] the 3 [lon,lat] geographic coordinates "      " [degrees]
            * xNbgh:  [(nT,3)   array of integers] the 3 IDs of the 3 neighbor triangles
        '''
        vpoints = np.unique( xPntID.flatten() )

        # Integers:
        self.nP        = len(vpoints) ; # number of points making the triangles
        self.nT        = nT ; # number of triangles (nT > nP)
        self.length    = nT ; #    "        "
        #
        # NumPy Arrays:
        self.PointIDs    = vpoints                                    ; # IDs of all the points making the triangles => shape = (nP)
        self.TriIDs      = np.array([i for i in range(nT)], dtype=int); # IDs of the triangles => shape = (nT)
        self.TriPointIDs = np.array(xPntID, dtype=int)                ; # 3 point IDs composing the triangles, CCW => shape = (nT,3)
        self.TriCoor     = np.array(xCoor)              ; # Coordinates of the 3 points composing the triangle => shape = (nT,3,2)
        self.neighbors   = np.array(xNbgh, dtype=int)
        #
        del vpoints

        
    def angles( self ):
        ''' Returns the (nT,3) array of the 3 angles (from 1st point to 3rd point, so counter-clockwize) '''
        return np.array( [ AnglesOfTriangle( self.TriCoor[i] )  for i in range(self.nT) ] )




class Quadrangle:

    def __init__( self, nQ, xPntID, xCoor ):
        '''
            * nQ:    [integer] number of Quadrangles
            * xPntID:  [(nQ,4)   array of integers] the 4 point IDs composing the quad, in counter-clockwize
            * xCoor: [(nQ,4,2) array of floats] the 4 [lon,lat] geographic coordinates "      " [degrees]
        '''
        vpoints = np.unique( xPntID.flatten() )
        
        # Integers:
        self.nP        = len(vpoints) ; # number of points making the quadrangles
        self.nQ        = nQ ; # number of quadrangles (nT > nP)        
        self.length    = nQ ; #    "        "

        # NumPy Arrays:
        self.PointIDs    = vpoints                                    ; # IDs of all the points making the quadrangles => shape = (nP)
        self.QuaIDs      = np.array([i for i in range(nQ)], dtype=int); # IDs of the quadrangles => shape = (nQ)
        self.QuaPointIDs = np.array(xPntID, dtype=int)                ; # 4 point IDs composing the quadrangles, CCW => shape = (nT,3)
        self.QuaCoor     = np.array(xCoor)              ; # Coordinates of the 3 points composing the quadrangle => shape = (nT,3,2)

        
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
            va1 = AnglesOfTriangle( np.array([ self.QuaCoor[i,j]  for j in [0,1,3] ]) ); # angles of triangle [ABD]
            va2 = AnglesOfTriangle( np.array([ self.QuaCoor[i,j]  for j in [2,3,1] ]) ); # angles of triangle [CDB]
            xA[i,:] = [ va1[0], va1[1]+va2[2], va2[0], va1[2]+va2[1] ]
            #
        return xA

