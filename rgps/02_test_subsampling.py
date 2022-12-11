#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
import numpy as np
from re import split

#from scipy.spatial import Delaunay

from gudhi import subsampling as sbspl


from climporn import epoch2clock
import mojito   as mjt

idebug=1

rzoom_fig = 5



if __name__ == '__main__':

    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <SELECTION_buoys_RGPS_streamXXX_XXX.npz>')
        exit(0)
    cf_npz = argv[1]

    #cfroot = str.replace( split('.npz',path.basename(cf_npz))[0] , 'SELECTION_buoys_RGPS_','' )
    #cf_npzT = './npz/T-mesh_'+cfroot+'.npz'
    #cf_npzQ = './npz/Q-mesh_'+cfroot+'.npz'
    #if (not path.exists(cf_npzT)) or (not path.exists(cf_npzQ)):
    # Have to build triangle and quadrangle mesh!
    #print('\n *** We are going to build triangle and quad meshes!')

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

    xCoor = np.array( [ [vx[i],vy[i]] for i in range(NbP) ] ) ; # original x,y cartesian cordinates of the RGPS data!


    # That's the original cloud of points, plotting it:

    kk = mjt.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='./figs/OriginalCloud.png',  lGeoCoor=False,
                         rangeX=[-1650,-700], rangeY=[-400,100], zoom=1.5 )
    
    # Projection, need to provide lon,lat, not distances:
    #kf = mjt.ShowBuoysMap( it, xCoor[:,0], xCoor[:,1], pvIDs=vIDs, cnmfig='OriginalCloud' )


    print(' *** Shape of xCoor: ',np.shape(xCoor))
    
    xx = sbspl.sparsify_point_set(xCoor, min_squared_dist=2500.)

    xx = np.array(xx)
    
    print(' *** Shape of xx: ',np.shape(xx))
    
    kk = mjt.ShowTQMesh( xx[:,0], xx[:,1], cfig='./figs/SparsifiedCloud.png',  lGeoCoor=False,
                         rangeX=[-1650,-700], rangeY=[-400,100], zoom=1.5 )


    
    #print(xx)
    
    #print(xCoor)
