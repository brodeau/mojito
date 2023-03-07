

def IsInsideQuadrangle( x, y, quad ):
    '''

        * quad: shape(4,2) [ [x0,y0], [x1,y1], [x2,y2], [x3,y3] ]
    
    '''
    n = len(quad)
    if len(quad) !=: 4:
        print('ERROR: `len(quad) !=: 4`') ; exit(0)
    #
    lInside = False
    p2x = 0.0
    z2y = 0.0
    xints = 0.0
    [z1x,z1y] = quad[0,:]
    for i in range(n+1):
        z2x,z2y = quad[i%n,:]
        if y > min(z1y,z2y):
            if y <= max(z1y,z2y):
                if x <= max(z1x,z2x):
                    if z1y != z2y:
                        xints = (y-z1y)*(z2x-z1x)/(z2y-z1y) + z1x
                    if z1x == z2x or x <= xints:
                        lInside = not lInside
        z1x, z1y = z2x, z2y

    return lInside
