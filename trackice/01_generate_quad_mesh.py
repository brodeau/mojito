#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#     CLIMPORN
#
#  Plot only the trajectories of all the buoys, even those who disapear
#
#  INPUT DATA: a `npz` file created with `trackice/scripts/traj2npz.py` (conversion from CSV to NPZ)
#
#    L. Brodeau, August 2022
#
# TO DO: use `nemo_box = cp.nemo_hbox(CNEMO,CBOX)` !!!
#
#
#  ABOUT input `npz` file:
#   * Name: should be of the form `NANUK4_ICE-BBM00_6h_19960101_19961031(_xxx).npz`
#
##################################################################

from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split

from scipy.spatial import Delaunay

import lbrgps   as lbr

from climporn import chck4f,Dates2NbDays

from netCDF4 import Dataset


idebug = 1

irec0 = 0 ; # record to start from...



narg = len(argv)
if not narg in [3,4]:
    print('Usage: '+argv[0]+' <file_trj.npz> <LSM_file> (iTsubsampl)')
    exit(0)

cf_npz = argv[1]
cf_lsm = argv[2]
# Subsampling in time...
itsubs = 1
if narg == 4 :
    itsubs = int(argv[3])

chck4f(cf_npz)
chck4f(cf_lsm)


#cnout  = str.replace( path.basename(cf_npz), '.npz', '.'+fig_type )


# Getting time info and time step from input npz file which is should look like NEMO output file:
vv = split('-|_', path.basename(cf_npz))

print(vv)


CCONF = vv[0]
print('\n *** CONF = '+CCONF)

cdt   = vv[3]
print('\n *** Time frequency = '+cdt)

cdt1, cdt2 = vv[4], vv[5]
print('\n *** Start and End dates => '+cdt1+' -- '+cdt2)

#cyear = cdt1[0:4]
#print('\n *** Year = '+cyear+'\n')

NbDays = Dates2NbDays(cdt1,cdt2)
print('     ==> number of days =', NbDays)

if cdt[-1]=='h':
    NbRecs = int(24/int(cdt[0:-1])*NbDays)
    if cdt=='1h': NbRecs = NbRecs+24 ; # fixme: why???? not fotr '6h' ???
else:
    print('ERROR: please adapt unit frequency: "'+cdt[-1]+'" !!!'); exit(0)
print('     ==> expected number of records =', NbRecs,'\n')


dir_conf = path.dirname(cf_npz)
if dir_conf == '':  dir_conf = '.'
print('\n *** dir_conf =',dir_conf,'\n')




#if not path.exists(cf_npz):

if not path.exists("npz"): mkdir("npz")


#############################3
print('\n *** Reading into '+cf_npz+' !!!')
with np.load(cf_npz) as data:
    NrTraj = data['NrTraj']
    xmask   = data['mask'][:,irec0]
    xIDs    = data['IDs'][:,irec0]
    xJIs    = data['JIs'][:,irec0]
    xJJs    = data['JJs'][:,irec0]
    #xFFs    = data['FFs'] ; # we do not care about the field...



if NrTraj != NbRecs-1:
    print('ERROR: NrTraj != NbRecs-1 !!!',NrTraj,NbRecs-1); exit(0)

print('\n *** Trajectories contain '+str(NrTraj)+' records...')
print('         => reading from record #'+str(irec0)+'!')

(NbBuoys,) = np.shape(xIDs)
print('\n *** There are '+str(NbBuoys)+' buoys at the begining...')


print('\n JIs =', xJIs[::10])
print('\n JIs =', xJJs[::10])

zlon, zlat = np.zeros(NbBuoys), np.zeros(NbBuoys)

# We need to load the NEMO's metric files to translate `jj,ji` to actual coordinates:
print('\n *** Reading "'+CCONF+'" metrics in "'+cf_lsm+'" ...')
with Dataset(cf_lsm) as id_lsm:
    #id_lsm = Dataset(cf_lsm)
    nb_dim = len(id_lsm.variables['tmask'].dimensions)
    Ni = id_lsm.dimensions['x'].size
    Nj = id_lsm.dimensions['y'].size
    print('    --- the shape of the '+CCONF+' domain appears to be Ni, Nj =', Ni, Nj)
    xlon_t = id_lsm.variables['glamt'][0,:,:]
    xlat_t = id_lsm.variables['gphit'][0,:,:]
    xlon_u = id_lsm.variables['glamu'][0,:,:]
    #xlat_u = id_lsm.variables['gphiu'][0,:,:]
    #xlon_v = id_lsm.variables['glamv'][0,:,:]
    xlat_v = id_lsm.variables['gphiv'][0,:,:]
    #if nb_dim == 4: XMSK  = id_lsm.variables[cnmsk][0,0,:,:] ; # t, y, x
    #if nb_dim == 3: XMSK  = id_lsm.variables[cnmsk][0,  :,:] ; # t, y, x
    #if nb_dim == 2: XMSK  = id_lsm.variables[cnmsk][    :,:] ; # t, y, x
    #if l_show_msh:
    #    Xlon = id_lsm.variables['glamu'][0,:,:]
    #    Xlat = id_lsm.variables['gphiv'][0,:,:]
    #id_lsm.close()
    print('      done.')


for jb in range(NbBuoys):
    rjj, rji = xJJs[jb], xJIs[jb]
    
    jj, rj = int(rjj)-1, rjj%1.   ; # F2C !
    ji, ri = int(rji)-1, rji%1.   ; # F2C !
    #print( ' jj, rj =', jj, rj)
    #print( ' ji, ri =', ji, ri)

    # #fixme: I'm still not sure whether 0.5 means T point or U,V point !!! => find out !!!

    if not rj in [0.,0.5] or not ri in [0.,0.5]:
        print('ERROR: not rj in [0.,0.5] or not ri in [0.,0.5] !!!')
        
    zlat[jb] = 2.*(0.5-rj)*xlat_t[jj,ji] + 2.*rj*xlat_v[jj,ji]
    zlon[jb] = 2.*(0.5-ri)*xlon_t[jj,ji] + 2.*ri*xlon_u[jj,ji]

    if idebug>0:
        print(' jb =',jb,': rjj, rji =',rjj, rji)
        print('      => zlat =', zlat[jb])
        print('      => zlon =', zlon[jb])

        
        print('')
#  -- frisrt working with geographic coordinates rather than cartesian coordinates...




exit(0)

# Generating triangular meshes out of the cloud of points:
TRI = Delaunay(xCoor)











exit(0)

cnmsk = 'tmask'
print('\n *** Reading "'+cnmsk+'" in meshmask file...')
id_lsm = Dataset(cf_lsm)
nb_dim = len(id_lsm.variables[cnmsk].dimensions)
Ni = id_lsm.dimensions['x'].size
Nj = id_lsm.dimensions['y'].size
#if i2 == 0: i2 = Ni
#if j2 == 0: j2 = Nj
if nb_dim == 4: XMSK  = id_lsm.variables[cnmsk][0,0,:,:] ; # t, y, x
if nb_dim == 3: XMSK  = id_lsm.variables[cnmsk][0,  :,:] ; # t, y, x
if nb_dim == 2: XMSK  = id_lsm.variables[cnmsk][    :,:] ; # t, y, x
if l_show_msh:
    Xlon = id_lsm.variables['glamu'][0,:,:]
    Xlat = id_lsm.variables['gphiv'][0,:,:]
id_lsm.close()
print('      done.')

print('\n The shape of the domain is Ni, Nj =', Ni, Nj)

