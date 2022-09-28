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

xTriangles = TRI.simplices.copy() ; # shape = (Nbt,3) A simplex of 2nd order is a triangle! *_*

(NbT,_) = nmp.shape(xTriangles) ; # NbT => number of triangles

xNeighbors = TRI.neighbors.copy() ;  # shape = (Nbt,3)

print('\n *** We have '+str(NbT)+' triangles!')

#print('\n',nmp.shape(TRI.points))
#exit(0)



if idebug>1:
    for jx in range(NbT):
        vpl = xTriangles[jx,:] ; # 3 point indices forming the triangle
        print(' Triangle #'+str(jx)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"')
        print('    => neighbor triangles are:',xNeighbors[jx,:],'\n')

if idebug>1:
    # In which simplex (aka triangle) is Grenoble:
    vlocate = nmp.array([(x_gre,y_gre)])
    kv_gre  = TRI.find_simplex(vlocate)
    jx      = kv_gre[0]
    vpl     = xTriangles[jx,:]
    print('\n *** Grenoble is located in triangle #'+str(jx)+':', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"\n')




# Show triangles on a map:
kk = lbr.ShowTMeshMap( Xcoor[:,0], Xcoor[:,1], xTriangles, cfig="01_Mesh_Map_TRIangles_Europe.png", pnames=vnam )


# Attempt to merge triangles into quadrangles:
# For now trying to merge #12 with #13
#   => aka "Andorra - Rome - Bern" with "Tunis - Rome - Andorra"
#  Each triangle inside the domain has 3 neighbors, so there are 3 options to merge
#  (provided the 3 neighbors are not already merged with someone!)


NbR = 0 ; # Number of quadrangles

idxTdead = []   ; # IDs of canceled triangles
idxTused  = []   ; # IDs of triangles already in use in a quadrangle
Quads       = []

for jT in range(NbT):
    # Loop along triangles

    v3pnts = xTriangles[jT,:] ; # 3 point IDs composing the triangle.

    if idebug>0:
        print('\n ******************************************************************')
        print(' *** Focus on triangle #'+str(jT)+' =>',[ vnam[i] for i in v3pnts ])
    
    if lbr.lTriangleOK(jT, xTriangles, Xcoor):        
        if idebug>0: print('       => disregarding this triangle!!! (an angle >120. or <30 degrees!)')
        idxTdead.append(jT) ; # Cancel this triangle
        #
    elif jT in idxTused:
        if idebug>0: print('       => this triangle is in use a quadrangle already defined!')
        #
        #
    else:
        # Triangle `jT` has a "decent" shape and has not been used to build a quad yet!
        # -----------------------------------------------------------------------------
        #
        vtmp   = xNeighbors[jT,:]
        vnghbs = vtmp[vtmp >= 0] ; # shrink it, only retain non `-1`-flagged values...
        NbN    = len(vnghbs)     ; # number of neighbors
        if idebug>0: print('       => its '+str(NbN)+' neighbor triangles are:', vnghbs)

        NgbrTvalid = [] ; # ID the valid neighbor triangles, i.e.: not dead, not already in use, and decent shape!
        for jN in vnghbs:
            lTok = (not jN in idxTdead)and(not jN in idxTused)and(not lbr.lTriangleOK(jN, xTriangles, Xcoor))
            if lTok:
                NgbrTvalid.append(jN)
                if idebug>1: print('          ==> triangle '+str(jN)+' is valid!')
            else:
                if idebug>1: print('          ==> triangle '+str(jN)+' is NOT valid!')
                    
        if len(NgbrTvalid)>0:
            # `jT` is a valid+available triangle with at least one valid neighbor in `NgbrTvalid`
            #   => need to check which of the neighbors in `NgbrTvalid` gives the best quadrangle!
            if idebug>0: print('       => valid neighbors for triangle #'+str(jT)+':',NgbrTvalid)

            #NNTok  = 0  ; # number of neighbor triangles are ok for forming a "decent" quad
            vjNok  = [] ; # stores the neighbor triangles that are ok for forming a "decent" quad
            vscore = [] ; #   => store their score of "OK-ness" !
            xidx   = []
            for jN in NgbrTvalid:
                if idebug>1: print('          ==> trying neighbor triangle '+str(jN)+':')
                vidx, vang  = lbr.QuadAnglesFrom2Tri( xTriangles, xNeighbors, jT, jN, Xcoor) #, pnam=vnam )
                lQok, score = lbr.lQuadOK( vang[:] )
                cc = 'does NOT'
                if lQok:
                    vjNok.append(jN)
                    vscore.append(score)
                    xidx.append(vidx)
                    cc = 'does'
                if idebug>1: print('            ===> "triangles '+str(jT)+'+'+str(jN)+'" '+cc+' give a valid Quad!')
            # Now we have to chose the best neighbor triangle to use (based on the score):
            if len(vjNok)>0:
                if idebug>0: print('       => We have '+str(len(vjNok))+' Quad candidates!')
                xidx = nmp.array(xidx)
                iwin = nmp.argmax(vscore)
                jN   = vjNok[iwin] ; # our winner neighbor triangle
                Quads.append(xidx[iwin,:]) ; # saving the winner Quad!
                idxTused.append(jT)
                idxTused.append(jN)
                if idebug>0: print('         ==> Selected Quad: "triangles '+str(jT)+'+'+str(jN)+'" give Quad',[vnam[i] for i in xidx[iwin,:]],'\n')
                #if len(vjNok)>1: exit(0) ; # LOLO DEBUG
        else:
            if idebug>0: print('       => No valid neighbors for this triangle...')


xQuads = nmp.array(Quads)
del Quads

(NbQ,_) = nmp.shape(xQuads)

if len(idxTused)/2 != NbQ or len(idxTused)%2 !=0:
    print('ERROR of agreement between number of merged triangles and quadrangles created!'); exit(0)

print('\n *** Triangles that have sucessfully be merged into acceptable Quads:\n   ==>', idxTused)

print('   ==> Summary of '+str(NbQ)+' generated quadrangles:')
for jQ in range(NbQ):
    print('    * Quad #'+str(jQ)+' => ', xQuads[jQ,:], '(', [ vnam[i] for i in xQuads[jQ,:] ],')')

# Show quadrangles on a map:
kk = lbr.ShowQMeshMap( Xcoor[:,0], Xcoor[:,1], xQuads, cfig="02_Mesh_Map_Quadrangles_Europe.png", pnames=vnam, TriMesh=xTriangles )


exit(0)











