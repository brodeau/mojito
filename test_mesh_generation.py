#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
#from os import path
import numpy as nmp

from scipy.spatial import Delaunay

import lbrgps   as lbr

idebug=2

y_gre, x_gre = 45.184369, 5.734251


xcities = [] ; ip = 0
#
xcities.append({"ID":ip,"city":"Paris"  , "lat":48.835334, "lon":2.353824 }); ip=ip+1
xcities.append({"ID":ip,"city":"Rome"   , "lat":41.89,     "lon":12.49    }); ip=ip+1
xcities.append({"ID":ip,"city":"Andorra", "lat":42.506939, "lon":1.521247 }); ip=ip+1
xcities.append({"ID":ip,"city":"Athen"  , "lat":37.984149, "lon":23.727984}); ip=ip+1
xcities.append({"ID":ip,"city":"Belgrad", "lat":44.817813, "lon":20.456897}); ip=ip+1
xcities.append({"ID":ip,"city":"Berlin" , "lat":52.517037, "lon":13.388860}); ip=ip+1
xcities.append({"ID":ip,"city":"Bern"   , "lat":46.948271, "lon":7.451451 }); ip=ip+1
xcities.append({"ID":ip,"city":"London" , "lat":51.510433, "lon":-0.129711}); ip=ip+1
xcities.append({"ID":ip,"city":"Esbjerg", "lat":55.477434, "lon":8.468160 }); ip=ip+1
xcities.append({"ID":ip,"city":"Brest",   "lat":48.389657, "lon":-4.481700}); ip=ip+1
xcities.append({"ID":ip,"city":"Tunis",   "lat":36.802481, "lon":10.168440}); ip=ip+1
xcities.append({"ID":ip,"city":"Madrid",  "lat":40.414060, "lon":-3.699336}); ip=ip+1
xcities.append({"ID":ip,"city":"Alger",   "lat":36.732591, "lon": 3.101878}); ip=ip+1
xcities.append({"ID":ip,"city":"Dublin",  "lat":53.341825, "lon":-6.267302}); ip=ip+1
xcities.append({"ID":ip,"city":"Porto",   "lat":41.144559, "lon":-8.622812}); ip=ip+1
xcities.append({"ID":ip,"city":"Tanger",  "lat":35.771494, "lon":-5.836354}); ip=ip+1
xcities.append({"ID":ip,"city":"Bergen",  "lat":60.378183, "lon": 5.333694}); ip=ip+1
xcities.append({"ID":ip,"city":"Stockholm","lat":59.310631,"lon":18.067388}); ip=ip+1
xcities.append({"ID":ip,"city":"Ghardaia","lat":32.483352,"lon":3.681163}); ip=ip+1

NbP = len(xcities) ; # number of points

print('\n *** We have '+str(NbP)+' cities!')


# Conversion of dictionary to Numpy arrays:
vIDs  = nmp.array([ xcities[jc]["ID"]                      for jc in range(NbP) ], dtype=int  )
Xcoor = nmp.array([[xcities[jc]["lon"],xcities[jc]["lat"]] for jc in range(NbP) ]             )
vnam  = nmp.array([ xcities[jc]["city"]                    for jc in range(NbP) ], dtype='U32')

if idebug>0:
    for jc in range(NbP):
        print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', lat='+str(round(Xcoor[jc,1],2))+', lon='+str(round(Xcoor[jc,0],2)))
    print('')


# Generating triangular meshes out of the cloud of points:
TRI = Delaunay(Xcoor)

Xtriangles = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

(NbT,_) = nmp.shape(Xtriangles) ; # NbT => number of triangles

Xneighbors = TRI.neighbors.copy() ;  # shape = (Nbt,3)

print('\n *** We have '+str(NbT)+' triangles!')

#print('\n',nmp.shape(TRI.points))
#exit(0)



if idebug>1:
    for jx in range(NbT):
        vpl = Xtriangles[jx,:] ; # 3 point indices forming the triangle
        print(' Triangle #'+str(jx)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"')
        print('    => neighbor triangles are:',Xneighbors[jx,:],'\n')

if idebug>1:
    # In which simplex (aka triangle) is Grenoble:
    vlocate = nmp.array([(x_gre,y_gre)])
    kv_gre  = TRI.find_simplex(vlocate)
    jx      = kv_gre[0]
    vpl     = Xtriangles[jx,:]
    print('\n *** Grenoble is located in triangle #'+str(jx)+':', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"\n')




# Show triangles on a map:

kk = lbr.ShowTMeshMap( Xcoor[:,0], Xcoor[:,1], Xtriangles, cfig="Mesh_Map_TRIangles_Europe.png", pnames=vnam )


# Attempt to merge triangles into quadrangles:
# For now trying to merge #12 with #13
#   => aka "Andorra - Rome - Bern" with "Tunis - Rome - Andorra"
#  Each triangle inside the domain has 3 neighbors, so there are 3 options to merge
#  (provided the 3 neighbors are not already merged with someone!)


NbR = 0 ; # Number of quadrangles

idxT_cancel = [] ; # IDs of canceled triangles


for jT in range(NbT):

    v3pnts = Xtriangles[jT,:] ; # 3 point IDs composing the triangle.
    #zcoorT = nmp.array([ Xcoor[j,:] for j in v3pnts ])    ; # Coordinates of the 3 points of a given triangle:
    if idebug>0: print('\n *** Focus on triange #'+str(jT)+' =>',[ vnam[i] for i in v3pnts ])
    
    if lbr.AnglesOfTriangleNotOK(jT, Xtriangles, Xcoor):
        # Cancel this triangle
        if idebug>0: print('       => disregarding this triangle!!! (an angle >120. or <30 degrees!)')
        idxT_cancel.append(jT)
        
    else:
        # Triangle seems fine!
        vtmp   = Xneighbors[jT,:]
        vnghbs = vtmp[vtmp >= 0] ; # shrink if `-1` are presents!
        NbN    = len(vnghbs)
        if idebug>0: print('       => its '+str(NbN)+' neighbor triangles are:', vnghbs)

        valid_nb = []
        for jN in vnghbs:
            if not jN in idxT_cancel:
                if idebug>1: print('          ==> triangle '+str(jN)+':',end='')
                if not lbr.AnglesOfTriangleNotOK(jN, Xtriangles, Xcoor):
                    valid_nb.append(jN)
                    if idebug>1: print(' is valid!')
                    
        if len(valid_nb)>0:
            # We have a valid (jT) triangle with at least one valid neighbor in `valid_nb`
            #  => need to check which of the neighbors in `valid_nb` gives the best quadrangle!
            if idebug>0: print('       => valid neighbors for triangle #'+str(jT)+':',valid_nb)
            for jN in valid_nb:
                if idebug>1: print('          ==> trying neighbor triangle '+str(jN)+':')
                vidx, vang = lbr.QuadAnglesFrom2Tri( Xtriangles, Xneighbors, jT, jN, Xcoor, pnam=vnam )
                print(vang)
                exit(0)

                print('')
        else:
            print('       => No valid neighbors for this triangle...')
        exit(0)
        # 1st neighbor:
        #Quad1 = Quad = lbr.Triangles2Quadrangle( Xtriangles, Xneighbors, jT, vnghbs[0], Xcoor, pnam=vnam )
exit(0)

#j1 = 12 ; j2 = 13 ; # Join 2 triangles with common segment: "Andorra-Rome"
j1 = 2  ; j2 = 5  ;  # Join 2 triangles with common segment: "London-Bergen"

Quad = lbr.Triangles2Quadrangle( Xtriangles, Xneighbors, j1, j2, Xcoor, pnam=vnam )

print('\n *** Quadran created by merging triangles '+str(j1)+' & '+str(j2)+':', '=>',[ vnam[i]         for i in Quad ] )
print('                             =>' ,[ (round(Xcoor[i,1],2),round(Xcoor[i,0],2)) for i in Quad ] )

