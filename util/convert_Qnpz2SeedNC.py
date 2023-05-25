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
iplot  = 1

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


    

    
    #cdate = QUADin.date
    

    
    ztim = QUADin.PointTime
    zIDs = QUADin.PointIDs
    zPXY = QUADin.PointXY

    idate = np.mean(ztim)
    cdate = mjt.epoch2clock(idate)

    print('   =>',Np,'points and',Nq,'Quads!\n      * Mean date: '+cdate)

    # Geographic coordinates:
    zPGC = mjt.CartNPSkm2Geo1D( zPXY, convArray='F' )

    if iplot>0:
        kf = mjt.ShowBuoysMap( idate, zPGC[:,0], zPGC[:,1], pvIDs=[], cfig='buoys_Q2NC.png', nmproj='CentralArctic', cnmfig=None,
                               ms=5, ralpha=0.5, lShowDate=True, zoom=1., title=None )



    
    ##############################
    # Time to save the stuff !!! #
    ##############################

    cf_nc_out = str.replace( cfQin, '.npz', '.nc' )

    print('   * Saving file: '+cf_nc_out)
    #kk = mjt.ncSaveCloudBuoys( cf_nc_out, ztim, vIDs, xYkm, xXkm, xlat, xlon, mask=xmsk,
    #                                   xtime=xtim, fillVal=mjt.FillValue, corigin='RGPS' )



    

    # Masking:
    #xXkm = np.ma.masked_where( xmsk==0, xXkm )
    #xYkm = np.ma.masked_where( xmsk==0, xYkm )
    
    #cout_root = 'SELECTION_RGPS_S'+'%3.3i'%(jS)+'_dt'+cdtbin_h


    # GENERATION OF COMPREHENSIVE NETCDF FILE:
    #  * 1 file per batch
    #  * for the time variable inside netCDF, we chose the mean time accross buoys in the bin used aka `vtim`
    #  * for the file name we chose the time bin bounds of scanned period
    #cdt1, cdt2 = split(':',mjt.epoch2clock(VT[0,1]))[0] , split(':',mjt.epoch2clock(VT[-1,2]))[0]  ; # keeps at the hour precision...
    #cdt1, cdt2 = str.replace( cdt1, '-', '') , str.replace( cdt2, '-', '')
    #cdt1, cdt2 = str.replace( cdt1, '_', 'h') , str.replace( cdt2, '_', 'h')




    #cf_nc_out = './nc/'+cout_root+'_'+cdt1+'_'+cdt2+'.nc'

