#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
#from os import path
import numpy as np

from climporn import epoch2clock
import mojito   as mjt

idebug=1

rzoom_fig = 3


if __name__ == '__main__':

    if not len(argv) in [3]:
        print('Usage: '+argv[0]+' <SELECTION_buoys_RGPS_streamXXX_XXX.npz> <min_dist_km>')
        exit(0)
    cf_npz = argv[1]
    cd_min = argv[2] ; rd_min = float(cd_min) ; cd_min = '%2.2i'%int(cd_min)

    print('\n *** Reading into '+cf_npz+' !!!')

    with np.load(cf_npz) as data:
        it     = data['itime']
        cdate  = str( data['date'] )
        Nbuoys = data['Npoints']
        vids   = data['vids']
        vx = data['vx']
        vy = data['vy']
        if len(vids) != len(vx) or len(vids) != len(vy):
            print('ERROR: len(vids) != len(vx) or len(vids) != len(vy)!') ; exit(0)

    NbP = len(vids) ; # number of points
    if NbP != Nbuoys: print('ERROR: NbP != Nbuoys !'); exit(0)

    if len(vx)!=NbP or len(vy)!=NbP:      print('ERROR Y13!') ; exit(0)
    if len(vids) != len(np.unique(vids)): print('ERROR Y14!') ; exit(0)

    print('\n *** Stream at '+epoch2clock(it)+' => '+str(NbP)+' points!')

    vIDs  = np.array( vids )
    del vids

    xCoor = np.array( [ [vx[i],vy[i]] for i in range(NbP) ] ) ; # original x,y cartesian cordinates of the RGPS data!


    # That's the original cloud of points, plotting it:

    kk = mjt.ShowTQMesh( xCoor[:,0], xCoor[:,1], cfig='./figs/01_OriginalCloud.png', ppntIDs=vIDs, lGeoCoor=False,
                         rangeX=[-1650,-700], rangeY=[-400,100], zoom=rzoom_fig )
    print(' *** Shape of xCoor: ',np.shape(xCoor))

    ########################################################################
    Nbuoys_ss, XY, vIDs_ss = mjt.SubSampCloud( rd_min, xCoor, vIDs )
    ########################################################################

    print(' *** New number of buoys =', Nbuoys_ss, ('(there was '+str(Nbuoys)+')'))
    
    kk = mjt.ShowTQMesh( XY[:,0], XY[:,1], cfig='./figs/02_SparsifiedCloud_'+cd_min+'km.png', ppntIDs=vIDs_ss,
                         lGeoCoor=False,
                         rangeX=[-1650,-700], rangeY=[-400,100], zoom=rzoom_fig )

