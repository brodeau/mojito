#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

#from sys import argv, exit
#from os import path
#import numpy as nmp

import pandas as pd


xcapitals = []
#
xcapitals.append({"city":"Paris"  , "lat":48.835334, "lon":2.353824 })
xcapitals.append({"city":"Rome"   , "lat":41.89,     "lon":12.49 })
xcapitals.append({"city":"Andorra", "lat":42.506939, "lon":1.521247 })
xcapitals.append({"city":"Athen"  , "lat":37.984149, "lon":23.727984})
xcapitals.append({"city":"Belgrad", "lat":44.817813, "lon":20.456897})
xcapitals.append({"city":"Berlin" , "lat":52.517037, "lon":13.388860})
xcapitals.append({"city":"Bern"   , "lat":46.948271, "lon":7.451451 })


Nbc = len(xcapitals)

print('\n We have '+str(Nbc)+' cities!')

print( xcapitals[2]["city"] )

