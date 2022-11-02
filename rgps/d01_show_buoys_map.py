#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
import numpy as np
import mojito as mjt


if not len(argv) in [3]:
    print('Usage: '+argv[0]+' <selection_file.npz> <zoom>')
    exit(0)
cf_npz = argv[1]
izoom  = int(argv[2])

print('')

cf_fig = str.replace( path.basename(cf_npz), '.npz', '.png' )

# Reading the quad meshes in the npz files:
with np.load(cf_npz) as data:
    itime = data['itime']
    vlon  = data['vlon']
    vlat  = data['vlat']
    vids  = data['vids']
    

# Show
mjt.ShowBuoysMap( itime, vlon, vlat, pvIDs=[], cfig=cf_fig, ms=5, ralpha=0.5 )
