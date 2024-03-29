#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
import mojito   as mjt

idebug=1

if not len(argv) in [3]:
    print('Usage: '+argv[0]+' <file_Q_mesh.npz> <zoom>')
    exit(0)
cf_npzQ = argv[1]
izoom   = int(argv[2])

print('')

cf_fig = str.replace( path.basename(cf_npzQ), '.npz', '.png' )

# Reading the quad meshes in the npz files:
QUA = mjt.LoadClassPolygon( cf_npzQ, ctype='Q' )

# Show the quads with only the points that define them:
kk = mjt.ShowTQMesh( QUA.PointXY[:,0], QUA.PointXY[:,1], cfig=cf_fig, QuadMesh=QUA.MeshVrtcPntIdx, lGeoCoor=False, zoom=izoom )
