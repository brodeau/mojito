import numpy as np
from itertools import combinations
from scipy.spatial import KDTree, ConvexHull, Delaunay


def sub_sample4(XY, d):
    """
    Sub-samples an array of (x,y) Cartesian coordinates so that points are at least d distance apart.

    Parameters:
        XY (ndarray): An array of shape (N, 2) containing the (x, y) coordinates of N points.
        d (float): The minimum distance between points in the sub-sampled array.

    Returns:
        ndarray: An array of shape (M, 2) containing the (x, y) coordinates of M points in the sub-sampled array.
    """
    # Convert the input array to a numpy array for convenience
    XY = np.asarray(XY)
    
    # Initialize the output array with the first point from the input array
    subsampled = [XY[0]]
    idxKeep = []    
    # Iterate over the input array, adding points that are at least d distance apart
    for i in range(1, len(XY)):
        # Check if the current point is at least d distance away from all points in the output array
        if np.all(np.linalg.norm(XY[i] - np.asarray(subsampled), axis=1) >= d):
            subsampled.append(XY[i])
            idxKeep.append(i)
    
    # Convert the output list to a numpy array and return it
    #    return np.asarray(subsampled)
    (nss,) = np.shape(idxKeep)
    # Convert the output list to a numpy array and return it
    #return np.asarray(subsampled)
    return nss, np.asarray(subsampled), np.asarray(idxKeep, dtype=int)




def get_quadrangular_mesh(XY):
    tri = Delaunay(XY)
    simplices = tri.simplices
    quadrangles = []
    for i in range(len(simplices)):
        neighbors = tri.neighbors[i]
        for j in range(3):
            if neighbors[j] != -1 and neighbors[j] > i:
                a = simplices[i, j]
                b = simplices[i, (j+1)%3]
                c = simplices[i, (j+2)%3]
                if tri.find_simplex(XY[a:b+1]) == i and tri.find_simplex(XY[b:c+1]) == i:
                    quadrangles.append([XY[a], XY[b], XY[c], XY[simplices[neighbors[j], :][tri.neighbors[neighbors[j]] != i][0]]])
    return np.array(quadrangles)


def quadrangle_area(p1, p2, p3, p4):
    '''Compute the area of a quadrangle defined by four points in the Cartesian plane.

    Parameters:
        p1, p2, p3, p4 (tuple): The (x,y) coordinates of the four points.

    Returns:
        float: The area of the quadrangle.
    '''
    # Compute the vectors of the diagonals
    d1 = np.array(p3) - np.array(p1)
    d2 = np.array(p4) - np.array(p2)
    
    # Compute the cross product of the vectors
    cp = np.cross(d1, d2)
    
    # Compute the magnitude of the cross product
    return abs(cp) / 2




def get_quadrangles(XY, Aref):
    '''Find quadrangles with an area close to a reference area in an array of x,y cartesian coordinates.

    Parameters:
        XY (numpy.ndarray): An array of shape (N, 2) containing the x,y cartesian coordinates of N points.
        Aref (float): The reference area for the quadrangles.

    Returns:
        numpy.ndarray: An array of shape (M, 8) containing the x,y cartesian coordinates of the M quadrangles,
        where each quadrangle is specified by the (x,y) coordinates of its four vertices in counter-clockwise order.
    '''
    # Compute the squared distances between all pairs of points
    D = np.sum((XY[:, np.newaxis, :] - XY[np.newaxis, :, :]) ** 2, axis=2)
    
    # Compute all combinations of four points
    combos = combinations(range(len(XY)), 4)
    
    # Initialize an empty list to store the quadrangles
    quads = []
    
    # Iterate over all combinations of four points
    for c in combos:
        # Check if the four points form a convex quadrangle
        p1, p2, p3, p4 = [XY[i] for i in c]
        if not is_convex(p1, p2, p3, p4):
            continue
        
        # Compute the area of the quadrangle
        #area = compute_area(p1, p2, p3, p4)
        area = quadrangle_area(p1, p2, p3, p4)
        
        # Check if the area is close to the reference area
        if not is_close(area, Aref):
            continue
        
        # Add the quadrangle to the list
        quads.append([p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1]])
    
    # Convert the list to a numpy array and return it
    return np.array(quads)

def is_convex(p1, p2, p3, p4):
    '''Check if four points form a convex quadrangle.'''
    # Compute the cross product of the vectors p1p2 and p1p3
    v1 = p2 - p1
    v2 = p3 - p1
    cp1 = v1[0] * v2[1] - v1[1] * v2[0]
    
    # Compute the cross product of the vectors p1p2 and p1p4
    v3 = p4 - p1
    cp2 = v1[0] * v3[1] - v1[1] * v3[0]
    
    # Compute the cross product of the vectors p3p4 and p3p2
    v4 = p2 - p3
    v5 = p4 - p3
    cp3 = v4[0] * v5[1] - v4[1] * v5[0]
    
    # Compute the cross product of the vectors p3p4 and p3p1
    v6 = p1 - p3
    cp4 = v4[0] * v6[1] - v4[1] * v6[0]
    
    # Check if the cross products have the same sign
    return (cp1 * cp2 > 0) and (cp3 * cp4 > 0)



#def compute_area(p1, p2, p3, p4):
#    '''Compute the area of
