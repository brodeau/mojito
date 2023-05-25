#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

'''
Takes a the npz file of quadrangles and convert it to a netCDF file fit to be used as the
seeding file for sitrack!
'''

from sys import argv, exit
from os import path, environ, makedirs
import numpy as np
from re import split

import mojito as mjt
from mojito import config as cfg

idebug = 0
iplot  = 1 ; lShowIDs=True ; zoom=5

#================================================================================================



if __name__ == '__main__':

    kk = cfg.initialize()

    for cd in ['npz','nc','figs']:
        makedirs( cd, exist_ok=True )
    #    if not path.exists('./'+cd): mkdir('./'+cd)
    #if not path.exists('./figs/SELECTION'): mkdir('./figs/SELECTION')

    ####################################################################################################
    narg = len(argv)
    if not narg in [2]:
        print('Usage: '+argv[0]+' <file_Quad.npz>)')
        exit(0)
    cfQin    = argv[1]
    #quality_mode = argv[5]
    
    #if not quality_mode  in ['thorough','model','rgps','xlose']:
    #    print('ERROR => unknow mode: '+mode+'!') ; exit(0)

    #idtbin_h = int(cdtbin_h)
    #

    mjt.chck4f(cfQin)

    QUADin = mjt.LoadClassPolygon( cfQin, ctype='Q' )


    Np = QUADin.nP
    Nq = QUADin.nQ


    

    
    
    ztim = QUADin.PointTime
    zIDs = QUADin.PointIDs
    zPXY = QUADin.PointXY

    idate = np.mean(ztim)
    cdate = mjt.epoch2clock(idate)

    print('   =>',Np,'points and',Nq,'Quads!\n      * Mean date: '+cdate)


    if idebug>0 and iplot>0:
        # Show the quads as just read:
        vrngX = mjt.roundAxisRange( QUADin.PointXY[:,0], rndKM=50. )
        vrngY = mjt.roundAxisRange( QUADin.PointXY[:,1], rndKM=50. )    
        kk = mjt.ShowTQMesh( QUADin.PointXY[:,0], QUADin.PointXY[:,1], cfig='QUADin_as_read.png',
                             ppntIDs=QUADin.PointIDs, QuadMesh=QUADin.MeshVrtcPntIdx, qIDs=QUADin.QuadIDs,
                             lGeoCoor=False, zoom=zoom, rangeX=vrngX, rangeY=vrngY, lShowIDs=lShowIDs )

    
    # Geographic coordinates:
    zPGC = mjt.CartNPSkm2Geo1D( zPXY, convArray='F' )

    if iplot>0:
        kf = mjt.ShowBuoysMap( idate, zPGC[:,0], zPGC[:,1], pvIDs=[], cfig='buoys_Q2NC.png', nmproj='CentralArctic', cnmfig=None,
                               ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )

    
    ##############################
    # Time to save the stuff !!! #
    ##############################

    cf_nc_out = str.replace( cfQin, 'npz', 'nc' )

    print('   * Saving file: '+cf_nc_out+' (only 1 record!)')
    kk = mjt.ncSaveCloudBuoys( cf_nc_out, [ idate ], zIDs, zPXY[:,1], zPXY[:,0], zPGC[:,1], zPGC[:,0],
                                       xtime=ztim, fillVal=mjt.FillValue, corigin='RGPS' )



    if idebug>0:
        # Control what we plan to do later on once the model has advected this cloud of points...
        
        # Now! We have a cloud of points read (saved) in the netCDF file that has all the points involved in the definition
        # of the quads we need.
        # We want to use the this cloud of point, together with the quads pattern known in the Quad file we use to reconstruct
        # exactly the same quads!!!
    
        xPxy2, vTime2, xQpnts2, vPids2, vQnam2, vQIDs2  = mjt.RecycleQuads( zPXY[:,:], ztim, zIDs, QUADin )
    
    
        kk = mjt.ShowTQMesh( xPxy2[:,0], xPxy2[:,1], cfig='QUAD_RECONSTRUCTED.png',
                             ppntIDs=vPids2, QuadMesh=xQpnts2, qIDs=vQIDs2,
                             lGeoCoor=False, zoom=zoom, rangeX=vrngX, rangeY=vrngY, lShowIDs=lShowIDs )
    
        #ppntIDs=vPids2, QuadMesh=QUADin.MeshVrtcPntIdx, qIDs=QUADin.QuadIDs,
    
