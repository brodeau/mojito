#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path
from glob import glob
import numpy as np
from re import split
from math import ceil,floor

#from scipy.spatial import Delaunay
import matplotlib.pyplot as plt


from climporn import epoch2clock, clock2epoch
import mojito   as mjt

idebug=1

cprefixIn='DEFORMATIONS_z' ; # Prefix of deformation files...

# Conversion from s-1 to day-1:
rconv = 24.*3600.


# Bin widths for pdfs
wbin_div = 0.0025 ; # day^-1
#wbin_div = 0.01 ; # day^-1



if not len(argv) in [2]:
    print('Usage: '+argv[0]+' <directory_input_npz_files>')
    exit(0)
cd_in = argv[1]
#cf_Q2 = argv[2]
#mrkrsz= int(argv[3])

# Polpulating deformation files available:
listnpz = np.sort( glob(cd_in+'/'+cprefixIn+'*.npz') )
nbFiles = len(listnpz)
print('\n *** We found '+str(nbFiles)+' deformation files into '+cd_in+' !')


kStreamName = np.zeros(nbFiles, dtype='U4')
kiDate      = np.zeros(nbFiles, dtype=int ) ; # date in epoch time
kNbPoints   = np.zeros(nbFiles, dtype=int ) ; # number of points in file

#kDate1      = np.zeros(nbFiles, dtype='U10')
kf = 0
for ff in listnpz:
    print('\n  # File: '+ff)
    fb = path.basename(ff)
    vf = split('_',fb)
    kStreamName[kf] = vf[2]
    #
    with np.load(ff) as data:
        idate = int( data['idate'] )
        #cdate = str( data['cdate'] )
        nPnts =      data['Npoints']
        #Xc    =      data['Xc']
        #Yc    =      data['Yc']
        #zdiv  =      data['divergence']
        #zshr  =      data['shear']

    kiDate[kf] = idate
    kNbPoints[kf] = nPnts

    print('   * Date = ',epoch2clock(kiDate[kf]))
    print('   * Nb. of points = ',kNbPoints[kf] )
        
    kf = kf+1

print('\n')

print('  ==> list of streams:', kStreamName[:])

nP = np.sum(kNbPoints)
print('  ==> Total number of points:', nP)



# Now that we know the total number of points we can allocate and fill arrays for divergence and shear
Zdiv = np.zeros(nP)
Zshr = np.zeros(nP)

jP = 0 ; # Counter from 0 to nP-1
kf = 0 ; # Counter for files, 0 to nbFiles-1
for ff in listnpz:
    jPe = jP + kNbPoints[kf]
    with np.load(ff) as data:
        zdiv  =      data['divergence']
        zshr  =      data['shear']
    #
    Zdiv[jP:jPe] = rconv*zdiv ; # day^-1
    Zshr[jP:jPe] = rconv*zshr ; # day^-1
    #
    jP = jPe
    kf = kf+1


#print(Zdiv)

div_min, div_max = np.min(Zdiv), np.max(Zdiv)
print(' min and max for div:', div_min, div_max)

vdmax =  ( round( 1.25 * max(abs(div_min),abs(div_max)), 2 ) )
print('    ==> x-axis max =',vdmax,' day^-1')

nBins = 2.*vdmax / wbin_div
if not nBins%1.==0.:
    print('ERROR: nBins is not an integer! nBins =',nBins)
    exit(0)

nBins = int(nBins)
print('nBins =',nBins)

xbin_bounds_div = [ -vdmax + float(i)*wbin_div for i in range(nBins+1) ]
xbin_bounds_div = np.round( xbin_bounds_div, 6 )
print('xbin_bounds_div =',xbin_bounds_div)

xbin_center_div = [ -vdmax+0.5*wbin_div + float(i)*wbin_div for i in range(nBins) ]
xbin_center_div = np.round( xbin_center_div, 6 )
print('xbin_center_div =',xbin_center_div)


PDF_div = np.zeros(nBins)

for iP in range(nP):
    rdiv = Zdiv[iP]

    #print('rdiv =',rdiv)

    jf = np.argmin( np.abs( xbin_center_div - rdiv ) )    
    #(idx1,) = np.where( xbin_bounds_div >=  rdiv )
    #(idx2,) = np.where( xbin_bounds_div <= rdiv )

    if not ( rdiv>=xbin_bounds_div[jf] and rdiv<xbin_bounds_div[jf+1] ):
        print(' Binning error!'); exit(0)

    PDF_div[jf] = PDF_div[jf]+1
    #print('jf = ', jf)
    #print('idx1 =',idx1)
    #print('idx2 =',idx2)

    #exit(0)



PDF_div[:] = PDF_div[:]/float(nP)




kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, name='Divergence', cfig='PDF_divergence.png', nx_subsamp=4  )


#kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, name='Divergence', cfig='PDF_divergence.png'  )



exit(0)
