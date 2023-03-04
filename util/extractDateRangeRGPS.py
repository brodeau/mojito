#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
###########################################################################

from sys import argv, exit
from os import path, environ, mkdir
import numpy as np

from re import split
import mojito as mjt

idebug = 0

cdt_pattern = 'YYYY-MM-DD_hh:mm:00' ; # pattern for dates

fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

list_expected_var = [ 'index', 'x', 'y', 'lon', 'lat', 'q_flag', 'time', 'idstream', 'streams' ]

FillValue = -9999.

#================================================================================================


if __name__ == '__main__':

    narg = len(argv)
    if not narg in [4]:
        print('Usage: '+argv[0]+' <file_RGPS.nc> <YYYYMMDD(_HH:MM)1> <YYYYMMDD(_HH:MM)2>')
        exit(0)
    cf_in    =     argv[1]
    cdate1   =     argv[2]
    cdate2   =     argv[3]
    ####################################################################################################
    ldp1, lhp1 = len(cdate1)==8, len(cdate1)==14
    ldp2, lhp2 = len(cdate2)==8, len(cdate2)==14

    cY1,  cY2  = cdate1[0:4], cdate2[0:4]
    cmm1, cmm2 = cdate1[4:6], cdate2[4:6]
    cdd1, cdd2 = cdate1[6:8], cdate2[6:8]

    chhmm1, chhmm2 = '00:00', '00:00'
    if lhp1: chhmm1 = cdate1[9:11]+':'+cdate1[12:14]
    if lhp2: chhmm2 = cdate2[9:11]+':'+cdate2[12:14]

    cdt1 = str.replace(cdt_pattern,'YYYY-MM-DD',cY1+'-'+cmm1+'-'+cdd1)
    cdt2 = str.replace(cdt_pattern,'YYYY-MM-DD',cY2+'-'+cmm2+'-'+cdd2)
    
    cdt1 = str.replace(cdt1,'hh:mm',chhmm1)
    cdt2 = str.replace(cdt2,'hh:mm',chhmm2)

    # File to save work in:
    cf_out = 'RGPS_'+split('_',cdt1)[0]+'_'+split('_',cdt2)[0]+'.nc4'
    
    print('\n *** Date range to restrain data to:')
    print(' ==> '+cdt1+' to '+cdt2 )

    idt1, idt2 = mjt.clock2epoch(cdt1), mjt.clock2epoch(cdt2)
    print( '   ===> in epoch time: ', idt1, 'to', idt2 )
    print( '       ====> double check: ', mjt.epoch2clock(idt1), 'to',  mjt.epoch2clock(idt2),'\n')

    # nP, nS, ztime, zy, zx, zlat, zlon, kBIDs, kStrm
    Np0, Ns0, ztime0, zykm0, zxkm0, zlat0, zlon0, zIDs0, zStrm0 = mjt.LoadDataRGPS( cf_in, list_expected_var )


    (idxKeep,) = np.where( (ztime0>=idt1) & (ztime0<=idt2) )
    
    
    kk = mjt.SaveRGPStoNC( cf_out, ztime0[idxKeep], zIDs0[idxKeep], zykm0[idxKeep], zxkm0[idxKeep],
                           zlat0[idxKeep], zlon0[idxKeep], Nstrm=Ns0, pSid=zStrm0[idxKeep] )
