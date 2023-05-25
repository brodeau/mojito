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
from glob import glob

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
        print('Usage: '+argv[0]+' <file_Quad_root*>)')
        exit(0)

    
    cQroot = argv[1]
    listnpz = np.sort( glob(cQroot+'*.npz') )

    Nrec = len(listnpz)

    print(' *** We have '+str(Nrec)+' files, and hence '+str(Nrec)+' records in output file!')
    for ff in listnpz:
        print('       * '+ff)
        mjt.chck4f(ff)


    for jr in range(Nrec):

        cfQin    = listnpz[jr]
        print('\n *** Record #'+str(jr+1)+':\n     => '+cfQin)

        # Dates as string from file name:
        cfQinB = path.basename(cfQin)
        print(cfQinB)
        zsfn = split('_',cfQinB)
                
        if jr==0:
            croot  = zsfn[0]+'_'+zsfn[1]+'_'+zsfn[2]+'_'+zsfn[3]
            cdt1   = str.replace(zsfn[-3],'-','h')[0:11]
            creskm = split('\.',zsfn[-1])[0]
        if jr==1:
            cdt2   = str.replace(zsfn[-2],'-','h')[0:11]
            
        
        QUADin = mjt.LoadClassPolygon( cfQin, ctype='Q' )

        Np = QUADin.nP
        Nq = QUADin.nQ

        ztim = QUADin.PointTime        
        zPXY = QUADin.PointXY

        if jr==0:
            zIDs = QUADin.PointIDs
            #
            vtim = np.zeros(    Nrec , dtype=int)
            xtim = np.zeros((Nrec,Np), dtype=int)
            xPXY = np.zeros((Nrec,Np,2))
            xPGC = np.zeros((Nrec,Np,2))
        
        idate = np.mean(ztim)
        cdate = mjt.epoch2clock(idate)        
        
        print('   =>',Np,'points and',Nq,'Quads!\n      * Mean date: '+cdate)

        if idebug>0 and iplot>0:
            # Show the quads as just read:
            vrngX = mjt.roundAxisRange( QUADin.PointXY[:,0], rndKM=50. )
            vrngY = mjt.roundAxisRange( QUADin.PointXY[:,1], rndKM=50. )    
            kk = mjt.ShowTQMesh( QUADin.PointXY[:,0], QUADin.PointXY[:,1], cfig='QUADin_as_read_rec%2.2i'%(jr)+'.png',
                                 ppntIDs=QUADin.PointIDs, QuadMesh=QUADin.MeshVrtcPntIdx, qIDs=QUADin.QuadIDs,
                                 lGeoCoor=False, zoom=zoom, rangeX=vrngX, rangeY=vrngY, lShowIDs=lShowIDs )
    
        
        # Geographic coordinates:
        zPGC = mjt.CartNPSkm2Geo1D( zPXY, convArray='F' )
        #print('shape(zPGC) =', np.shape(zPGC))
        #exit(0)
        
        if iplot>0:
            kf = mjt.ShowBuoysMap( idate, zPGC[:,0], zPGC[:,1], pvIDs=[], cfig='buoys_Q2NC_rec%2.2i'%(jr)+'.png',
                                   nmproj='CentralArctic', cnmfig=None,
                                   ms=5, ralpha=0.5, lShowDate=True, zoom=1, title=None )

        #kk = mjt.ncSaveCloudBuoys( cf_nc_out, [ idate ], zIDs, zPXY[:,1], zPXY[:,0], zPGC[:,1], zPGC[:,0],
        #                                   xtime=ztim, fillVal=mjt.FillValue, corigin='RGPS' )

            

        # Storing for current record:
        vtim[jr]     = idate
        xPXY[jr,:,:] = zPXY[:,:]
        xPGC[jr,:,:] = zPGC[:,:]
        xtim[jr,:]   = ztim[:]
            

        if idebug>0:

            # Control what we plan to do later on once the model has advected this cloud of points...
            
            # Now! We have a cloud of points read (saved) in the netCDF file that has all the points involved in the definition
            # of the quads we need.
            # We want to use the this cloud of point, together with the quads pattern known in the Quad file we use to reconstruct
            # exactly the same quads!!!
        
            xPxy2, vTime2, xQpnts2, vPids2, vQnam2, vQIDs2  = mjt.RecycleQuads( xPXY[jr,:,:], xtim[jr,:], zIDs, QUADin )
        
            kk = mjt.ShowTQMesh( xPxy2[:,0], xPxy2[:,1], cfig='QUAD_RECONSTRUCTED_rec%2.2i'%(jr)+'.png',
                                 ppntIDs=vPids2, QuadMesh=xQpnts2, qIDs=vQIDs2,
                                 lGeoCoor=False, zoom=zoom, rangeX=vrngX, rangeY=vrngY, lShowIDs=lShowIDs )



    del QUADin


        
    ##############################
    # Time to save the stuff !!! #
    ##############################

    #cf_nc_out = str.replace( cfQin, 'npz', 'nc' )


    # File must end: ${YEAR}????h??_${YEAR}????h??${cxtraRES}${XTRASFX}.nc

    #print('cdt1 =',cdt1)
    #print('cdt2 =',cdt2)    
    #print('creskm =',creskm)

    cf_nc_out = './nc/'+croot+'_'+cdt1+'_'+cdt2+'_'+creskm+'.nc'
    
    
    print('   * Saving file: '+cf_nc_out+' ('+str(Nrec)+' record!)')

    kk = mjt.ncSaveCloudBuoys( cf_nc_out, vtim, zIDs, xPXY[:,:,1], xPXY[:,:,0], xPGC[:,:,1], xPGC[:,:,0],
                               xtime=xtim, fillVal=mjt.FillValue, corigin='RGPS' )
    



            
    
        
            #ppntIDs=vPids2, QuadMesh=QUADin.MeshVrtcPntIdx, qIDs=QUADin.QuadIDs,
    
