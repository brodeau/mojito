#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
#from os import path, mkdir
import numpy as np

#from mojito import initialize, updateConfig4Scale
#import mojito as cfg
#import config as cfg
from mojito import config as cfg

kk = cfg.initialize()

print(' rc_day2sec =', cfg.rc_day2sec)


for ks in [ 10, 20, 40, 80, 160, 320, 640 ]:

    kk = cfg.updateConfig4Scale( ks, 6 )
