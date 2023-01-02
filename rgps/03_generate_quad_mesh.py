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

idebug=2

# Selection of appropriate quadrangles:
rTang_min =  10. ; # minimum angle tolerable in a triangle [degree]
rTang_max = 120. ; # maximum angle tolerable in a triangle [degree]
#
#rQang_min =  65.  ; # minimum angle tolerable in a quadrangle [degree]
#rQang_max = 115.  ; # maximum angle tolerable in a quadrangle [degree]
rQang_min =  60.  ; # minimum angle tolerable in a quadrangle [degree]
rQang_max = 120.  ; # maximum angle tolerable in a quadrangle [degree]
rdRatio_max = 0.7 ; # value that `1 - abs(L/H)` should not overshoot!

rL_nom     = 10. ; # nominal length [km] of the sides of the quadrangle we intend to build

rzoom_fig = 5



if __name__ == '__main__':

    if not len(argv) in [2,3]:
        print('Usage: '+argv[0]+' <SELECTION_buoys_RGPS_streamXXX_XXX.npz> (<min_pt_spacing_km>)')
        exit(0)
    cf_npz = argv[1]    
    l_force_min_scale = ( len(argv) == 3 and argv[2] != '0' )
    


    rftol = 0.3 ; #fixme
    
    if l_force_min_scale:
        cd_min = argv[2]
        rd_min = float(cd_min)
        #rftol  = 0.5 ; #fixme !!! should rather take into account the deviation from 10km...
        rftol = 0.3 * rd_min/rL_nom
        #
        rL_nom = rd_min
        


    cL_nom = '%2.2i'%int(round(rL_nom,0))
        
    # Strings for names of output files:
    cfroot = str.replace( split('.npz',path.basename(cf_npz))[0] , 'SELECTION_buoys_RGPS_','' )
    if l_force_min_scale:
        cfroot += '_'+cL_nom+'km_Sampled'
    else:
        cfroot += '_'+cL_nom+'km_NoSample'
    
    # #fixme: move the 2 coeffs to header...
    rf1 , rf2 = 1.-rftol , 1.+rftol
    rQarea_min = rf1*rL_nom*rL_nom  ; # min area allowed for Quadrangle [km^2]
    rQarea_max = rf2*rL_nom*rL_nom  ; # max area allowed for Quadrangle [km^2]
    
    cf_npzT = './npz/T-mesh_'+cfroot+'.npz'
    cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'

    if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ)):

        print('\n *** We are going to build triangle and quad meshes!')
        print('     => desired scale for quadrangles = '+str(int(rL_nom))+'km')
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

        itt     = 0
        l_happy = False
        rfact_corr = 1.
        rdev = 1.
        
        while not l_happy:
            itt = itt + 1
            
            # Just prior to Delaunay we may have to sub-sample in space the cloud of point
            if l_force_min_scale:
    
                if idebug>1:
                    kk = mjt.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='./figs/00_Original_'+cfroot+'.png',
                                         pnames=vPnam, ppntIDs=vIDs,
                                         lGeoCoor=False, zoom=rzoom_fig )
                                         #lGeoCoor=False, rangeX=[-1650,-700], rangeY=[-400,100], zoom=rzoom_fig )            
                # Update!!!! #fixme
                rd_min = rfact_corr * rL_nom

                print('\n *** Applying spatial sub-sampling with radius rd_min = '+str(round(rd_min,2))+' triangles!')
                
                NbPss, zCoor, zIDs, zPnam = mjt.SubSampCloud( rd_min, xCoor, vIDs,  pNames=vPnam ) ; #lilo
                
                if idebug>1:
                    kk = mjt.ShowTQMesh( zCoor[:,0], zCoor[:,1], cfig='./figs/00_SubSamp_'+cfroot+'.png',
                                         pnames=vPnam, ppntIDs=vIDs,
                                         lGeoCoor=False, zoom=rzoom_fig )
                                         #lGeoCoor=False, rangeX=[-1650,-700], rangeY=[-400,100], zoom=rzoom_fig )
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
            if len(xQpnts)<=0: exit(0)
    
            (NbQ,_) = np.shape(xQpnts)
            print('\n *** We have '+str(NbQ)+' quadrangles!')
    
            # Conversion to the `Quadrangle` class (+ we change IDs from triangle world [0:nT] to that of quad world [0:nQ]):
            QUADS = mjt.Quadrangle( xQcoor, xQpnts, vPQids, vQnam, date=cdate )

            del xQpnts, xQcoor, zCoor
            
            zsides = QUADS.lengths()
            rl_average_side = np.mean(zsides)
            if idebug>0:

                zareas = QUADS.area()
                rl_average_area = np.mean(zareas)
                print('    ==> average side length is '+str(round(rl_average_side,3))+' km')
                print('    ==> average area is '+str(round(rl_average_area,1))+' km^2')
                del zareas, zsides

            rdev_old = rdev
            rdev = rl_average_side - rL_nom
            l_happy = ( abs(rdev) < 0.05 ) ; # average quadrangle side is close to expected nominal scale

            if not l_happy:
                # Linear fit of actual correction as a function of `rL_nom`
                rfc = 0.008*rL_nom + 0.56                
                if itt==1: ralpha = (1.-rfc) / rdev ; # equivalent to a correction of `rfc`
                if itt>1 and copysign(1,rdev) == -copysign(1,rdev_old):
                    ralpha = ralpha/2. ; # change of sign of deviation => we decrease alpha!
                # will go for a next round with a correction factor becoming increasingly smaller than 1:
                rfact_corr = min( max(0.6 , rfact_corr - ralpha * rdev ) , 0.95 )
                print(' +++++ NEW ralpha, rfact_corr = ',ralpha, rfact_corr)
                
        ###################
        # #while not l_happy
        
        # Save the triangular mesh info:
        mjt.SaveClassPolygon( cf_npzT, TRIAS, ctype='T' )

        # Save the quadrangular mesh info:
        mjt.SaveClassPolygon( cf_npzQ, QUADS, ctype='Q' )

        del TRIAS, QUADS

    #if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ))
    ############################################################


    # Reading the triangle and quad class objects in the npz files:
    TRI = mjt.LoadClassPolygon( cf_npzT, ctype='T' )
    QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )


    if not path.exists('./figs'): mkdir('./figs')

    # Show triangles on a map:
    kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/01_Mesh_Triangles_'+cfroot+'.png',
                         TriMesh=TRI.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig)

    # Show triangles together with the quadrangles on a map:
    kk = mjt.ShowTQMesh( TRI.PointXY[:,0], TRI.PointXY[:,1], cfig='./figs/02_Mesh_Quadrangles_'+cfroot+'.png',
                         TriMesh=TRI.MeshVrtcPntIdx,
                         pX_Q=QUA.PointXY[:,0], pY_Q=QUA.PointXY[:,1], QuadMesh=QUA.MeshVrtcPntIdx,
                         lGeoCoor=False, zoom=rzoom_fig)

    ## Show only points composing the quadrangles:
    #kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/03_Mesh_Points4Quadrangles_'+cfroot+'.png',
    #                     lGeoCoor=False )

    # Show only the quads with only the points that define them:
    kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig='./figs/03_Mesh_Points4Quadrangles_'+cfroot+'.png',
                         QuadMesh=QUA.MeshVrtcPntIdx, lGeoCoor=False, zoom=rzoom_fig)


    print('\n *** rfact_corr was:',rfact_corr)
