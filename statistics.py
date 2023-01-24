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

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

# Conversion from s-1 to day-1:
rconv = 24.*3600.


# Bin widths for pdfs
#wbin_div = 0.0025 ; # day^-1
#wbin_div = 0.001 ; # day^-1
#wbin_div = 0.0005 ; # day^-1
wbin_div = 0.00025 ; # day^-1
#wbin_div = 0.0001 ; # day^-1
max_div = 0.5  ; # day^-1
xdiv_rng=[-0.04,0.04] ; # x-range we want on the x-axis of the plot


#wbin_shr = 0.0025 ; # day^-1
#wbin_shr = 0.001 ; # day^-1
#wbin_shr = 0.0005 ; # day^-1
wbin_shr = 0.00025 ; # day^-1
#wbin_shr = 0.0001 ; # day^-1
max_shr = 1. ; # day^-1
xshr_rng=[0.,0.1] ; # x-range we want on the x-axis of the plot

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
kiDate      = np.zeros(nbFiles, dtype=int ) ; # date in epoch time at which deformations were calculated
kNbPoints   = np.zeros(nbFiles, dtype=int ) ; # number of points in file

#kDate1      = np.zeros(nbFiles, dtype='U10')
kf = 0
for ff in listnpz:
    print('\n  # File: '+ff)
    fb = path.basename(ff)
    vf = split('_',fb)
    kStreamName[kf] = vf[1]
    #
    with np.load(ff) as data:
        rdate = int( data['time'] )
        #cdate = str( data['cdate'] )
        nPnts =      data['Npoints']
        #Xc    =      data['Xc']
        #Yc    =      data['Yc']
        #zdiv  =      data['divergence']
        #zshr  =      data['shear']

    kiDate[kf] = rdate
    kNbPoints[kf] = nPnts

    print('   * Stream: '+kStreamName[kf] )
    print('   * Date = ',epoch2clock(kiDate[kf]))
    print('   * Nb. of points = ',kNbPoints[kf] )
        
    kf = kf+1

print('\n')

#print('  ==> list of streams:', kStreamName[:])

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

if not max_div:
    div_min, div_max = np.min(Zdiv), np.max(Zdiv)
    print(' min & max for div:', div_min, div_max)
    max_div =  ( round( 1.25 * max(abs(div_min),abs(div_max)), 2 ) )
print('    ==> x-axis max =',max_div,' day^-1')

nBinsD = 2.*max_div / wbin_div
if not nBinsD%1.==0.:
    print('ERROR: nBinsD is not an integer! nBinsD =',nBinsD); exit(0)
nBinsD = int(nBinsD)
#print('nBinsD =',nBinsD)

xbin_bounds_div = [ -max_div + float(i)*wbin_div for i in range(nBinsD+1) ]
xbin_bounds_div = np.round( xbin_bounds_div, 6 )
#print('xbin_bounds_div =',xbin_bounds_div)

xbin_center_div = [ -max_div+0.5*wbin_div + float(i)*wbin_div for i in range(nBinsD) ]
xbin_center_div = np.round( xbin_center_div, 6 )
#print('xbin_center_div =',xbin_center_div)



# For the shear:
if not max_shr:
    shr_max = np.max(Zshr)
    print(' max for shr:', shr_max)
    max_shr =  ( round(1.25*shr_max, 2) )
print('    ==> x-axis max =',max_shr,' day^-1')

nBinsS = max_shr / wbin_shr
if not nBinsS%1.==0.:
    print('ERROR: nBinsS is not an integer! nBinsS =',nBinsS); exit(0)
nBinsS = int(nBinsS)
#print('nBinsS =',nBinsS)

xbin_bounds_shr = [  float(i)*wbin_shr for i in range(nBinsS+1) ]
xbin_bounds_shr = np.round( xbin_bounds_shr, 6 )
#print('xbin_bounds_shr =',xbin_bounds_shr)

xbin_center_shr = [ 0.5*wbin_shr + float(i)*wbin_shr for i in range(nBinsS) ]
xbin_center_shr = np.round( xbin_center_shr, 6 )
#print('xbin_center_shr =',xbin_center_shr)

#exit(0)




PDF_div = np.zeros(nBinsD)
PDF_shr = np.zeros(nBinsS)

for iP in range(nP):
    rdiv = Zdiv[iP]
    jf = np.argmin( np.abs( xbin_center_div - rdiv ) )    
    if not ( rdiv>=xbin_bounds_div[jf] and rdiv<xbin_bounds_div[jf+1] ):
        print(' Binning error on divergence!')
        print('  => divergence =',rdiv)
        print('  => bounds =',xbin_bounds_div[jf],xbin_bounds_div[jf+1])
        exit(0)
    PDF_div[jf] = PDF_div[jf]+1

    rshr = Zshr[iP]
    jf = np.argmin( np.abs( xbin_center_shr - rshr ) )    
    if not ( rshr>=xbin_bounds_shr[jf] and rshr<xbin_bounds_shr[jf+1] ):
        print(' Binning error on shear!')
        print('  => shear =',rshr)
        print('  => bounds =',xbin_bounds_shr[jf],xbin_bounds_shr[jf+1])
        exit(0)
    PDF_shr[jf] = PDF_shr[jf]+1



PDF_div[:] = PDF_div[:]/float(nP)
PDF_shr[:] = PDF_shr[:]/float(nP)

kk = mjt.PlotPDFdef( xbin_bounds_div, xbin_center_div, PDF_div, Np=nP, name='Divergence', cfig='PDF_divergence.png',
                     xrng=xdiv_rng, wbin=wbin_div, title='RGPS' )

kk = mjt.PlotPDFdef( xbin_bounds_shr, xbin_center_shr, PDF_shr, Np=nP, name='Shear', cfig='PDF_shear.png',
                     xrng=xshr_rng, wbin=wbin_shr, title='RGPS' )

