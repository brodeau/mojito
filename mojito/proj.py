
import numpy as np


Pi = 3.141592653589793
EarthRad = 6300. ;

toRad = Pi/180.

def NorthStereoProj( pphi, plam, lam0=0., phi0=90. ):
    #
    # https://mathworld.wolfram.com/StereographicProjection.html
    #
    # to km (since Earth radius provided in km)!!!
    #
    
    zsinPhi0  = np.sin(toRad*phi0)
    zsinPhi   = np.sin(toRad*pphi)
    zcosPhi0  = np.cos(toRad*phi0)
    zcosPhi   = np.cos(toRad*pphi)    
    zcosLambL = np.cos(toRad*(plam - lam0))
    

    zk = 2.*EarthRad / ( 1. + zsinPhi0*zsinPhi + zcosPhi0*zcosPhi*zcosLambL )

    zx = zk * zcosPhi*np.sin(toRad*(plam - lam0))

    zy = zk * ( zcosPhi0*zsinPhi - zsinPhi0*zcosPhi*zcosLambL )


    return zy, zx

    






    
