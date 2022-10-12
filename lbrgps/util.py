#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import exit
#from os import path, mkdir
import numpy as np


def LoadDist2CoastNC( cNCfile ):
    from climporn import chck4f
    from netCDF4  import Dataset
    #
    print('\n *** [util.LoadDist2CoastNC()] Loading "distance to coast" from file:')
    print('      '+cNCfile)
    chck4f(cNCfile)
    with Dataset(cNCfile) as id_in:
        vlon  = id_in.variables['lon'][:]
        vlat  = id_in.variables['lat'][:]
        xdist = id_in.variables['dist'][:,:]
        print('       => ok!\n')
        #
    return vlon, vlat, xdist


def Dist2Coast( lon0, lat0, plon, plat, pdist2coat ):
    '''
       Returns the distance to the nearest coast of a given point (lon,lat)
        INPUT:
          * lon0, lat0: coordinates (scalars) [degrees East], [degrees West]
          * plon:       1D array of longitudes assosiated to `pdist2coat`
          * plat:       1D array of latitudes  assosiated to `pdist2coat`
          * pdist2coat: 2D array containing rasterized distance to coast (originally read in NC file) [km]
    '''
    from climporn import degE_to_degWE
    #
    rx = np.mod( lon0, 360. ) ; # [0:360] frame
    vx = np.mod( plon, 360. ) ; # [0:360] frame
    # Are we dangerously close to the [0--360] cut?
    #  => then will work in the [-180:180] frame:
    if rx > 355.:
        rx = degE_to_degWE( rx )
        vx = degE_to_degWE( vx )
        #print(' lon0, lat0 =', rx, lat0)
    ip = np.argmin( np.abs(  vx[:] - rx  ) )
    jp = np.argmin( np.abs(plat[:] - lat0) )
    #print(' ip, jp =', ip, jp)
    #print(' Nearest lon, lat =', plon[ip], plat[jp])
    del vx, rx
    return max( pdist2coat[jp,ip] , 0. )



def OrderCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''

    # sort the points based on their x-coordinates
    isortX  = np.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]

    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]

    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = np.argsort(leftMost[:,1])
    leftMost   = leftMost[isortL,:]
    (tl, bl)   = leftMost
    [i1l, i2l] = isortML[isortL]

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = np.argsort(rghtMost[:,1])
    rghtMost   = rghtMost[isortR,:]
    (tr, br)   = rghtMost
    [i1r, i2r] = isortMR[isortR]

    #idx = np.concatenate([idxL,idxR])
    #isort = np.array([isortX[i] for i in np.concatenate([isortL,isortR])])
    #print('LOLO, idx=', idx)

    del isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR

    # return the coordinates in top-left, top-right,
    # bottom-right, and bottom-left order
    return np.array([tl, tr, br, bl], dtype="float32") , np.array([i1l, i1r, i2r, i2l])


def OrderCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''
    zz = xcoor.copy()
    iz = np.array([0,1,2,3])
    zt, it = OrderCW(xcoor)
    #
    zz[1:4,:] = zt[:0:-1,:]
    iz[1:4]   = it[:0:-1]
    del zt, it
    return zz, iz



def SortIndicesCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter-clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''
    # sort the points based on their x-coordinates
    isortX  = np.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]
    #print('LOLO: Sorted by longitude => isortX =', isortX)

    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]

    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = np.argsort(leftMost[:,1])
    [i1l, i2l] = isortML[isortL]
    #print('LOLO: isortL, idxL =', isortL, [i1l, i2l] )

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = np.argsort(rghtMost[:,1])
    [i1r, i2r] = isortMR[isortR]
    #print('LOLO: isortR, idxR =',isortR, [i1r, i2r] )

    #idx = np.concatenate([idxL,idxR])
    #isort = np.array([isortX[i] for i in np.concatenate([isortL,isortR])])
    #print('LOLO, idx=', idx)

    del xSorted, isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR

    return np.array([i1l, i1r, i2r, i2l])


def mergeNPZ( list_npz_files, cf_out='merged_file.npz' ):
    '''
       Merge several npz files of the same date into a single one!
    '''

    
