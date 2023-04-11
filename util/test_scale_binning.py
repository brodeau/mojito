#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
#from os import path, mkdir
import numpy as np

from mojito import updateConfig4Scale




for ks in [ 10, 20, 40, 80, 160, 320, 640 ]:

    kk = updateConfig4Scale( ks, 6 )
