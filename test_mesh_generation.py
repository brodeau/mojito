#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

#from sys import argv, exit
#from os import path
import numpy as nmp

from scipy.spatial import Delaunay

#import pandas as pd

import cartopy.crs as ccrs
#import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#from mpl_toolkits.basemap import Basemap



xcapitals = []
#
xcapitals.append({"city":"Paris"  , "lat":48.835334, "lon":2.353824 })
xcapitals.append({"city":"Rome"   , "lat":41.89,     "lon":12.49 })
xcapitals.append({"city":"Andorra", "lat":42.506939, "lon":1.521247 })
xcapitals.append({"city":"Athen"  , "lat":37.984149, "lon":23.727984})
xcapitals.append({"city":"Belgrad", "lat":44.817813, "lon":20.456897})
xcapitals.append({"city":"Berlin" , "lat":52.517037, "lon":13.388860})
xcapitals.append({"city":"Bern"   , "lat":46.948271, "lon":7.451451 })
# 51.510433, -0.129711
xcapitals.append({"city":"London" , "lat":51.510433, "lon":-0.129711 })
# 55.477434, 8.468160
xcapitals.append({"city":"Esbjerg", "lat":55.477434, "lon":8.468160 })
# 48.389657, -4.481700
xcapitals.append({"city":"Brest",   "lat":48.389657, "lon":-4.481700 })
# 36.802481, 10.168440
xcapitals.append({"city":"Tunis",   "lat":36.802481, "lon":10.168440 })
# 40.414060, -3.699336
xcapitals.append({"city":"Madrid",  "lat":40.414060, "lon":-3.699336 })


Nbc = len(xcapitals)

print('\n We have '+str(Nbc)+' cities!')

print( xcapitals[0]["city"] )

yk, xk = nmp.zeros(Nbc), nmp.zeros(Nbc)

for jc in range(Nbc):
    yk[jc], xk[jc] =  xcapitals[jc]["lat"], xcapitals[jc]["lon"]

#print(yk) ; print('') ; print(xk)

X1  = nmp.vstack((xk,yk)).T
tri = Delaunay(X1)



#--- grafics -------
figX = 25; figY = 18
fig1 = plt.figure(figsize=(figX, figY), facecolor='white')

myProjection = ccrs.PlateCarree()
ax = plt.axes(projection=myProjection)
ax.stock_img()
ax.set_extent([-25, 40, 35, 65], crs=myProjection)

plt.triplot(X1[:,0], X1[:,1], tri.simplices.copy(), color='r', linestyle='-',lw=2)
plt.plot(X1[:,0], X1[:,1], 's', color='w')

plt.scatter(xk,yk,s=1000,c='w')
#for i, txt in enumerate(city):
#    ax.annotate(txt, (X1[i,0], X1[i,1]), color='k', fontweight='bold')

plt.savefig('Europe_A.png')
#plt.show()
