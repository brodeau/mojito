#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
#from os import path
import numpy as nmp

#from scipy.spatial import Delaunay

import lbrgps   as lbr


vangles_perfect = nmp.array([ 90. , 90., 90., 90. ])

vangles_good = nmp.array([ 89. , 92., 85., 93. ])

vangles_bad = nmp.array([ 122. , 105., 77., 82. ])


vangles_terrible = nmp.array([ 132. , 115., 62., 72. ])


lok, scr = lbr.lQuadOK( vangles_perfect )
print("For perfect we have:", lok, scr )


lok, scr = lbr.lQuadOK( vangles_good )
print("For good we have:", lok, scr )

lok, scr = lbr.lQuadOK( vangles_bad )
print("For bad we have:", lok, scr )

lok, scr = lbr.lQuadOK( vangles_terrible )
print("For terrible we have:", lok, scr )
