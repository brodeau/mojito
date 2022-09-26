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


xcapitals = [] ; ip = 0
#
xcapitals.append({"ID":ip,"city":"Paris"  , "lat":48.835334, "lon":2.353824 }); ip=ip+1
xcapitals.append({"ID":ip,"city":"Rome"   , "lat":41.89,     "lon":12.49    }); ip=ip+1
xcapitals.append({"ID":ip,"city":"Andorra", "lat":42.506939, "lon":1.521247 }); ip=ip+1
xcapitals.append({"ID":ip,"city":"Athen"  , "lat":37.984149, "lon":23.727984}); ip=ip+1
xcapitals.append({"ID":ip,"city":"Belgrad", "lat":44.817813, "lon":20.456897}); ip=ip+1
xcapitals.append({"ID":ip,"city":"Berlin" , "lat":52.517037, "lon":13.388860}); ip=ip+1
xcapitals.append({"ID":ip,"city":"Bern"   , "lat":46.948271, "lon":7.451451 }); ip=ip+1
xcapitals.append({"ID":ip,"city":"London" , "lat":51.510433, "lon":-0.129711}); ip=ip+1
xcapitals.append({"ID":ip,"city":"Esbjerg", "lat":55.477434, "lon":8.468160 }); ip=ip+1
xcapitals.append({"ID":ip,"city":"Brest",   "lat":48.389657, "lon":-4.481700}); ip=ip+1
xcapitals.append({"ID":ip,"city":"Tunis",   "lat":36.802481, "lon":10.168440}); ip=ip+1
xcapitals.append({"ID":ip,"city":"Madrid",  "lat":40.414060, "lon":-3.699336}); ip=ip+1
xcapitals.append({"ID":ip,"city":"Alger",   "lat":36.732591, "lon": 3.101878}); ip=ip+1
# 53.341825, -6.267302
xcapitals.append({"ID":ip,"city":"Dublin",  "lat":53.341825, "lon":-6.267302}); ip=ip+1
# 41.144559, -8.622812
xcapitals.append({"ID":ip,"city":"Porto",   "lat":41.144559, "lon":-8.622812}); ip=ip+1
# 35.771494, -5.836354
xcapitals.append({"ID":ip,"city":"Tanger",  "lat":35.771494, "lon":-5.836354}); ip=ip+1
# 60.378183, 5.333694
xcapitals.append({"ID":ip,"city":"Bergen",  "lat":60.378183, "lon": 5.333694}); ip=ip+1
# 59.310631, 18.067388
xcapitals.append({"ID":ip,"city":"Stockholm","lat":59.310631,"lon":18.067388}); ip=ip+1


Nbc = len(xcapitals)

print('\n *** We have '+str(Nbc)+' cities!')


# Conversion of dictionary to Numpy vectors:
vIDs =  nmp.zeros(Nbc, dtype=int)
vlat, vlon = nmp.zeros(Nbc), nmp.zeros(Nbc)
vnam = nmp.zeros(Nbc, dtype='U32')

for jc in range(Nbc):
    vIDs[jc] = int( xcapitals[jc]["ID"] )
    vlat[jc], vlon[jc] =  xcapitals[jc]["lat"], xcapitals[jc]["lon"]
    vnam[jc] = xcapitals[jc]["city"]

if idebug>0:
    for jc in range(Nbc):
        print(' * '+vnam[jc]+': ID='+str(vIDs[jc])+', lat='+str(round(vlat[jc],2))+', lon='+str(round(vlon[jc],2)))
    print('')
    
#print(vlat) ; print('') ; print(vlon)
#print(vIDs)
#exit(0)

X1  = nmp.vstack((vlon,vlat)).T ; # Concatenate `vlon` (1D) and `vlat` (1D) in a single 2D array and transpose

TRI = Delaunay(X1)


Xsimplices = TRI.simplices.copy() ; # A simplex of 2nd order is a triangle! *_*


(Nbt,np) = nmp.shape(Xsimplices) ; # Nbt => number of traingles | np has to be = 3 !!! (triangles!)

print('\n *** We have '+str(Nbt)+' triangles!')

#if len(Xsimplices[:,0]) != Nbc: print('ERROR Z0!'); exit(0)

if idebug>0:
    for jx in range(Nbt):
        vpl = Xsimplices[jx,:] ; # 3 point indices forming the triangle
        print(' Triangle #'+str(jx)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"')
        #print('     '+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]])
        print('    => neighbor triangles are:',TRI.neighbors[jx],'\n')


# In which simplex (aka triangle) is Grenoble:
vlocate = nmp.array([(x_gre,y_gre)])
kv_gre= TRI.find_simplex(vlocate)
jx    = kv_gre[0]
vpl = Xsimplices[jx,:]
print('\n *** Grenoble is located into triangle #'+str(jx)+': ', vpl[:],'aka "'+vnam[vpl[0]]+' - '+vnam[vpl[1]]+' - '+vnam[vpl[2]]+'"\n')




# Show triangles on a map:
kk = lbr.ShowMeshMap( X1[:,0], X1[:,1], Xsimplices, cfig="Mesh_Map_TRIangles_Europe.png", pnames=vnam )

