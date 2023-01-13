#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

# TO DO: add standard deviation on the analysis of quadrangles properties...


from sys import argv, exit
from os import path
import numpy as np
from re import split
from math import copysign
from scipy.spatial import Delaunay

from climporn import epoch2clock
import mojito   as mjt

idebug=1

# Selection of appropriate quadrangles:
rTang_min =  10. ; # minimum angle tolerable in a triangle [degree]
rTang_max = 120. ; # maximum angle tolerable in a triangle [degree]
#
rQang_min =  65.  ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 115.  ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.7 ; # value that `1 - abs(L/H)` should not overshoot!


rcAtol = 0.3333 ; # coefficient of tolerance for the acceptation of the area of the quadrangles

rtol        = 0.3 ; # +- tolerance in [km] to accept a given scale. Ex: 15.19 km is accepted for 15 km !!!
rd_nom_data = 10. ; # default/nominal point spacing in [km] of the data

rzoom_fig = 5



if __name__ == '__main__':

    if not len(argv) in [2,3]:
        print('Usage: '+argv[0]+' <SELECTION_buoys_RGPS_streamXXX_XXX.npz> (<min_pt_spacing_km>)')
        exit(0)
    cf_npz = argv[1]    
    l_force_min_scale = ( len(argv) == 3 and argv[2] != '0' )
    
    rd_spacing = rd_nom_data
    if l_force_min_scale:
        rd_spacing = float(argv[2])
        if rd_spacing < rd_nom_data:
            print('ERROR: you cannot use a spacing smaller than that of the data!'); exit(0)
    cL_spacing = '%2.2i'%int(round(rd_spacing,0))


    print('\n\n**********************************************************************')
    print(' *** Ok! Scale we are going to use is: `rd_spacing` = ',str(rd_spacing)+'km')
    
    # Strings for names of output files:
    cfroot = str.replace( split('.npz',path.basename(cf_npz))[0] , 'SELECTION_buoys_RGPS_','' )
    if l_force_min_scale:
        cfroot += '_'+cL_spacing+'km_Sampled'
    else:
        cfroot += '_'+cL_spacing+'km_NoSample'
    
    rA_nom = rd_spacing*rd_spacing ; # expected nominal area of the quadrangles [km^2]
    rf1 , rf2 = 1.-rcAtol , 1.+rcAtol
    rQarea_min = rf1*rA_nom  ; # min area allowed for Quadrangle [km^2]
    rQarea_max = rf2*rA_nom  ; # max area allowed for Quadrangle [km^2]

    print('    => tolerance on the QUADs is: +-'+str(round(rcAtol*100.,2))+'%')
    print('    ==> '+str(rA_nom)+'km^2 ('+str(rQarea_min)+'km^2 < A < '+str(rQarea_max)+'km^2)')
    print('**********************************************************************\n')

    
    l_someQuads = True
    cf_npzT = './npz/T-mesh_'+cfroot+'.npz'
    cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'

    if not path.exists(cf_npzT) or not path.exists(cf_npzQ) or idebug>0:

        print('\n *** We are going to build triangle and quad meshes!')
        print('     => desired scale for quadrangles = '+str(int(rd_spacing))+'km')
        print('     => area range for a quadrangle to qualify: '+str(int(rQarea_min))+'km^2 < A < '+str(int(rQarea_max))+'km^2')
        
        print('\n *** Reading into '+cf_npz+' !!!')
        with np.load(cf_npz) as data:
            it     = data['itime']
            cdate  = str( data['date'] )
            Nbuoys = data['Npoints']
            vids   = data['vids']
            vx = data['vx']
            vy = data['vy']
            if len(vids) != len(vx) or len(vids) != len(vy): print('ERROR Y11!') ; exit(0)

        NbP = len(vids) ; # number of points
        if NbP != Nbuoys: print('ERROR: NbP != Nbuoys !'); exit(0)

        if len(vx)!=NbP or len(vy)!=NbP:      print('ERROR Y13!') ; exit(0)
        if len(vids) != len(np.unique(vids)): print('ERROR Y14!') ; exit(0)

        print('\n *** Stream at '+epoch2clock(it)+' => '+str(NbP)+' points!')

        vIDs  = np.array( vids )
        del vids


        # Name for each point:
        vPnam = np.array( [ str(i) for i in vIDs ], dtype='U32' )


        xCoor = np.array( [ [vx[i],vy[i]] for i in range(NbP) ] ) ; # original x,y cartesian cordinates of the RGPS data!
        #                                                           # => Polar Stereographic proj., lon_0=-45, lat_ts=70        
        if idebug>2:
            for jc in range(NbP):
                print(' * Name: "'+vPnam[jc]+'": ID='+str(vIDs[jc])
                      +', x ='+str(round(xCoor[jc,0],2))+', y ='+str(round(xCoor[jc,1],2)))
            print('')


        
        l_happy = False
        itt     = 0
        rfcorr  = 1.
        rdev    = 1.
        
        while not l_happy:
            itt = itt + 1
            
            # Just prior to Delaunay we may have to sub-sample in space the cloud of point
            if l_force_min_scale:
                    
                rd_ss = rfcorr * rd_spacing ; # Correct `rd_spacing` to get closer to requested radius (based on QUADs to be generated)

                print('\n *** Applying spatial sub-sampling! Threshold radius: '+str(round(rd_ss,2))+'km')                
                NbPss, zCoor, zIDs, zPnam = mjt.SubSampCloud( rd_ss, xCoor, vIDs,  pNames=vPnam )
                
                if idebug>0:
                    # Shows the cloud of buoys (with buoys' IDs) on the Cartesian plane (km)
                    # After and before subsampling
                    zrx = [ np.min(xCoor[:,0])-50. , np.max(xCoor[:,0])+50. ]
                    zry = [ np.min(xCoor[:,1])-50. , np.max(xCoor[:,1])+50. ]
                    #
                    kk = mjt.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='./figs/00_Original_'+cfroot+'.png',
                                         ppntIDs=vIDs, lGeoCoor=False, zoom=rzoom_fig, rangeX=zrx, rangeY=zry )
                    kk = mjt.ShowTQMesh( zCoor[:,0], zCoor[:,1], cfig='./figs/00_SubSamp_'+cfroot+'.png',
                                         ppntIDs=zIDs, lGeoCoor=False, zoom=rzoom_fig, rangeX=zrx, rangeY=zry )
                #lilo
            else:
                NbPss = NbP
                zCoor = xCoor
                zIDs  = vIDs
                zPnam = vPnam
                
    
            # Generating triangular meshes out of the cloud of points:
            TRI = Delaunay(zCoor)
    
            xTpnts = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*
    
            (NbT,_) = np.shape(xTpnts) ; # NbT => number of triangles
    
            xNeighborIDs = TRI.neighbors.copy() ;  # shape = (Nbt,3)
    
            print('\n *** We have '+str(NbT)+' triangles!')
    
            # Conversion to the `Triangle` class:
            TRIAS = mjt.Triangle( zCoor, xTpnts, xNeighborIDs, zIDs, zPnam ) ; #lolo
    
            del xTpnts, xNeighborIDs, TRI
            
            
            # Merge triangles into quadrangles:
            xQcoor, vPQids, xQpnts, vQnam = mjt.Tri2Quad( TRIAS, iverbose=idebug, anglRtri=(rTang_min,rTang_max),
                                                          ratioD=rdRatio_max, anglR=(rQang_min,rQang_max),
                                                          areaR=(rQarea_min,rQarea_max) )
            l_someQuads = (len(xQpnts)>0)
            #if len(xQpnts)<=0: exit(0)

            if l_someQuads:
                (NbQ,_) = np.shape(xQpnts)
                print('\n *** We have '+str(NbQ)+' quadrangles!')
    
                # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
                QUADS = mjt.Quadrangle( xQcoor, xQpnts, vPQids, vQnam, date=cdate )

                del xQpnts, xQcoor, zCoor
            
                zsides = QUADS.lengths()
                zareas = QUADS.area()
                rl_average_side = np.mean(zsides)
                rl_average_scal = np.mean( np.sqrt(zareas) )
                rl_average_area = np.mean(zareas)
                print('    ==> average scale (sqrt[A]) is '+str(round(rl_average_scal,3))+' km')
                print('    ==> average side length is '+str(round(rl_average_side,3))+' km')
                print('    ==> average area is '+str(round(rl_average_area,1))+' km^2')
                del zareas, zsides
    
                if l_force_min_scale:
                    rdev_old = rdev
                    #rdev = rl_average_side - rd_spacing
                    rdev = rl_average_scal - rd_spacing
                    l_happy = ( abs(rdev) < rtol ) ; # average quadrangle side is close to expected nominal scale

                    if not l_happy and itt==5:
                        # We give up after 5 itterations!
                        print(' +++++ WE GIVE UP !!! ++++++')
                        l_happy = True
                    
                    if not l_happy:
                        # Linear fit of actual correction as a function of `rd_spacing`
                        rfc = 0.008*rd_spacing + 0.56                
                        if itt==1: ralpha = (1.-rfc) / rdev ; # equivalent to a correction of `rfc`
                        if itt>1 and copysign(1,rdev) == -copysign(1,rdev_old):
                            ralpha = ralpha/1.3 ; # change of sign of deviation => we decrease alpha!
                        # will go for a next round with a correction factor becoming increasingly smaller than 1:
                        #rfcorr = min( max(0.6 , rfcorr - ralpha * rdev ) , 0.95 )
                        rfcorr = min( max(0.5 , rfcorr - ralpha * rdev ) , 0.95 )
                        print(' +++++ NEW itteration with: ralpha, rfcorr = ',ralpha, rfcorr)
                    #
                else:
                    l_happy = True ; # we work with the data's nominal scale so we just go on here...
            else:
                l_happy = True ; # we coud not build any QUAD so we move on anyways...
                
        ###################
        # #while not l_happy

        # Save the triangular mesh info:
        mjt.SaveClassPolygon( cf_npzT, TRIAS, ctype='T' )
        del TRIAS

        if l_someQuads:
            # Save the quadrangular mesh info:
            mjt.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q' )
            del QUADS            

    #if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ))
    ############################################################


    # Reading the triangle and quad class objects in the npz files:
    TRI = mjt.LoadClassPolygon( cf_npzT, ctype='T' )
    if l_someQuads: QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )


    if not path.exists('./figs'): mkdir('./figs')

    # Show triangles on a map:
    kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/01_Mesh_Triangles_'+cfroot+'.png',
                         TriMesh=TRI.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig, rangeX=zrx, rangeY=zry)

    if l_someQuads:
        # Show triangles together with the quadrangles on a map:
        kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/02_Mesh_Quadrangles_'+cfroot+'.png',
                             TriMesh=TRI.MeshVrtcPntIdx,
                             pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshVrtcPntIdx,
                             lGeoCoor=False, zoom=rzoom_fig, rangeX=zrx, rangeY=zry)
    
        ## Show only points composing the quadrangles:
        #kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/03_Mesh_Points4Quadrangles_'+cfroot+'.png',
        #                     lGeoCoor=False )
    
        # Show only the quads with only the points that define them:
        kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/03_Mesh_Points4Quadrangles_'+cfroot+'.png',
                             QuadMesh=QUA.MeshVrtcPntIdx, ppntIDs=QUA.PointIDs, lGeoCoor=False, zoom=rzoom_fig, rangeX=zrx, rangeY=zry)


    print('\n *** `rfcorr` was:',rfcorr)

    if not l_someQuads:
        print('\n *** NO QUADs could be built!!!   :(\n')