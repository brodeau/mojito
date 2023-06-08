#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
from glob import glob
import numpy as np
from re import split
#
import mojito   as mjt
from mojito import config as cfg

idebug=1
iplot=1

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...


def lIntersect2Quads( pQxy1, pQxy2 ):
    '''
       * pQxy1: cartesian coordinates of a quad, shape = (4,2)
       * pQxy2: cartesian coordinates of another quad, shape = (4,2)

       => will return true if the 2 quads intersect each other
    '''
    from shapely.geometry import Polygon
    #
    (n4,n2) = np.shape(pQxy1)
    if n4!=4 or n2!=2:
        print('ERROR [lIntersect2Quads()]: wrong shape for `pQxy1` !'); exit(0)
    if np.shape(pQxy2) != (n4,n2):
        print('ERROR [lIntersect2Quads()]: wrong shape for `pQxy2` !'); exit(0)
    #
    zP1 = Polygon( [ (pQxy1[i,0],pQxy1[i,1]) for i in range(4) ] )
    zP2 = Polygon( [ (pQxy2[i,0],pQxy2[i,1]) for i in range(4) ] )
    #
    return zP1.intersects(zP2)



def lOverlapAnyOtherQuad( pQxy, pMQ,  iverbose=0 ):
    '''
       * pQxy: cartesian coordinates of 1 quad, shape = (4,2)
       * pMQ:  cartesian coordinates of a groupd of quads, shape = (nq,4,2)

       => will return true if the quad `pQxy` overlaps partially (or completely) at least 1 quad of `pMQ` !
    '''
    from shapely.geometry import Polygon
    #
    (n4,n2) = np.shape(pQxy)
    if n4!=4 or n2!=2:
        print('ERROR [lOverlapAnyOtherQuad()]: wrong shape for `pQxy` !'); exit(0)
    (nq,n4,n2) = np.shape(pMQ)
    if n4!=4 or n2!=2:
        print('ERROR [lOverlapAnyOtherQuad()]: wrong shape for `pMQ` !'); exit(0)
    if iverbose>0: print(' *** [lOverlapAnyOtherQuad()]: we have',nq,'quads to test against...')
    #
    li = False
    for jq in range(nq):
        li = lIntersect2Quads( pQxy, pMQ[jq,:,:] )
        if li:
            break
    #
    return li




if __name__ == '__main__':

    kk = cfg.initialize()

    if not len(argv) in [6]:
        print('Usage: '+argv[0]+' <directory_input_npz_files> <batch.N> <dtbin_h> <creskm> <string_id_origin>')
        exit(0)

    cd_in  = argv[1]
    kbatch = argv[2]
    cdtbin = argv[3]
    creskm = argv[4]
    cidorg = argv[5]

    sbatch = 'S%3.3i'%(int(kbatch))



    print(cd_in+'/'+cprefixIn+'*'+cidorg+'_'+sbatch+'_'+cdtbin+'*_'+creskm+'km.npz')
    listnpz1 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'_'+sbatch+'_'+cdtbin+'*_'+creskm+'km.npz') )
    listnpz2 = np.sort( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'_'+sbatch+'_'+cdtbin+'*_*-'+creskm+'km.npz') )

    if len(listnpz1)>0 and len(listnpz2)>0:
        print('ERROR: we have both npz files with suffixes lile `*_Xkm` and `_Y-Xkm` !!!'); exit(0)
    if len(listnpz2)>0:
        lrlstKM = True
        listnpz = listnpz2
    else:
        lrlstKM = FalsenbF
        listnpz = listnpz1


    # Polpulating deformation files available:
    nbF = len(listnpz)
    print('\n *** We found '+str(nbF)+' deformation files into '+cd_in+' !')
    print('      => files to gather into a single one:',listnpz,'\n')

    # First, populate the number of deformation quadrangles in each file:
    kNbPoints = np.zeros(nbF, dtype=int ) ; # number of points in file
    vsubRes   = np.zeros(nbF, dtype=int ) ; # coarsifying resolution


    kf = 0
    for cf in listnpz:

        css   = split( '-', split( '_', split('\.',cf)[0] )[-1] )[0]
        vsubRes[kf] = int(css)
        csskm = css+'km'

        with np.load(cf) as data:
            nPnts =      data['Npoints']
            if kf==0:
                corigin = data['origin']
                reskm = data['reskm_nmnl']
        print('\n  # File: '+cf+' ('+csskm+') => '+str(nPnts)+' deformation quads!')
        kNbPoints[kf] = nPnts
        kf+=1


    # Indices of files we keep a the end:
    idxF2K   = [] ; # indices of files we keep a the end.
    idxF_Q2K = [] ; # for each retained file, the indices of quads to retain

    # Winner:
    kw = np.argmax(kNbPoints)
    Nmax = kNbPoints[kw]
    (idxmax,) = np.where(kNbPoints==Nmax)
    nn = len(idxmax)
    # Takes the one the closest to the reference resolution
    if nn>1:
        ki = np.argmin( np.abs(vsubRes[idxmax] - int(creskm)) )
        print(vsubRes[idxmax] - (int(creskm)-2), ki) ; # "-2" to be closer to the actual scale used to coarsify!
        kw = idxmax[ki]
    Nmax = kNbPoints[kw]
    print(' *** Winner is at index:',kw,'(',Nmax,'points),',vsubRes[kw],'km,',nn,'files with this number of points.')
    #
    idxF2K.append(kw)
    idxQ2K = np.arange(Nmax, dtype=int)
    idxF_Q2K.append( idxQ2K )


    # The "in use" quads coordinates array:
    nQmaxF = nbF*Nmax ; # the ever max number of gathered quads we can expect, final version of `nQmaxF` will have less points...
    zInUseXY = np.zeros((nQmaxF,4,2)) - 9999.
    zmsk     = np.zeros( nQmaxF , dtype='i1' )

    # 1/ Save all the quads of winner:
    cf = listnpz[kw]
    print(' *** will keep all the quads found in file:',cf)
    with np.load(cf) as data:
        zInUseXY[0:Nmax,:,0], zInUseXY[0:Nmax,:,1] = data['X4'], data['Y4']
        zmsk[0:Nmax] = 1
    #
    knxt = Nmax
    print('')


    print(' Shape of zInUseXY =',np.shape(zInUseXY))


    kf = 0
    for cf in listnpz:

        if kf != kw:

            cf = listnpz[kf]
            print('\n *** Getting the quads found in file:',cf)
            with np.load(cf) as data:
                zX4, zY4 = data['X4'], data['Y4']
            #
            (nQ,_) = np.shape(zX4)
            print('     ==>',nQ,'quads...')

            icpt = 0
            idxQ2K = []
            for jQ in range(nQ):
                zQxy = np.array([zX4[jQ,:],zY4[jQ,:]]).T ; # shape = (4,2)  (1 quad => x,y coordinates for 4 points)
                if not lOverlapAnyOtherQuad( zQxy, zInUseXY[:knxt,:,:] ):
                    print('       * quad #',jQ,'is SELECTED! ('+str(vsubRes[kf])+'km)')
                    zInUseXY[knxt,:,:] = zQxy[:,:]
                    zmsk[knxt] = 1
                    if icpt==0:
                        idxF2K.append(kf) ; # at least 1 quad selected, we keep this file!
                    idxQ2K.append(jQ)
                    icpt += 1
                    knxt += 1
                else:
                    print('       * quad #',jQ,' did not make it...')

            idxQ2K = np.array(idxQ2K, dtype=int)
            (nQslctd,) = np.shape(idxQ2K)
            if nQslctd>0:
                # we found at least 1 quad to keep:
                print('     => '+str(nQslctd)+' quads selected!' )
                idxF_Q2K.append( idxQ2K )

        ### if kf != kw
        kf += 1

    ### for cf in listnpz

    idxF2K = np.array( idxF2K, dtype=int )

    nQtot = np.sum(zmsk)



    print('\n\n   ### SUMMARY ###\n    At the end we selected '+str(nQtot)+' quads!')
    print('     ==> retained files and indices of quads to keep:')
    jf = 0
    qcpt = 0
    for cf in listnpz[idxF2K]:
        print('    - '+cf+' ('+str(vsubRes[idxF2K[jf]])+'km) => indices of quads to keep:',idxF_Q2K[jf])
        qcpt += len(idxF_Q2K[jf])
        jf += 1
    print('')

    if qcpt != nQtot:
        print('ERROR: something went wrong! `qcpt != nQtot`'); exit(0)
    
    #zX4 = np.array(zX4)
    #print(np.shape(zX4))



    # Allocation for the arrays we will save in the new npz file
    zXc, zYc = np.zeros(nQtot)    , np.zeros(nQtot)
    zX4, zY4 = np.zeros((nQtot,4)), np.zeros((nQtot,4))
    zdiv, zshr, ztot = np.zeros(nQtot), np.zeros(nQtot), np.zeros(nQtot)
    zAq = np.zeros(nQtot)

    jf = 0
    jq = 0
    for cf in listnpz[idxF2K]:
        #print('    - '+cf+' ('+str(vsubRes[idxF2K[jf]])+'km) => indices of quads to keep:',idxF_Q2K[jf])
        ikeep = idxF_Q2K[jf]
        nQr = len(ikeep)
        jqe = jq+nQr
        #
        with np.load(cf) as data:
            zXc[jq:jqe] = data['Xc'][ikeep]
            zYc[jq:jqe] = data['Yc'][ikeep]
            zX4[jq:jqe,:] = data['X4'][ikeep,:]
            zY4[jq:jqe,:] = data['Y4'][ikeep,:]
            zdiv[jq:jqe] = data['divergence'][ikeep]
            zshr[jq:jqe] = data['shear'][ikeep]
            ztot[jq:jqe] = data['total'][ikeep]
            zAq[jq:jqe] = data['quadArea'][ikeep]
            #
            if jf==0:
                itimeC = data['time']
                ctimeC = data['date']
        #        
        jf += 1
        jq = jqe
    ### cf in listnpz[idxF2K]
    print('')


    # Save the deformation data:
    #cfout = './npz/DEFORMATIONS_'+cfnm+'.npz'
    cfout = 'TEST.npz'
    
    np.savez_compressed( cfout, time=itimeC, date=ctimeC, Npoints=nQtot,
                         Xc=zXc, Yc=zYc, X4=zX4, Y4=zY4, divergence=zdiv, shear=zshr, total=ztot,
                         quadArea=zAq, origin=corigin, reskm_nmnl=reskm )


    
    exit(0)

















    ######






    kBatchName  = np.zeros(nbF, dtype='U4')
    kiDate      = np.zeros(nbF, dtype=int ) ; # date in epoch time at which deformations were calculated
    kNbPoints   = np.zeros(nbF, dtype=int ) ; # number of points in file

    list_date = []
    kf = 0
    for cf in listnpz:
        print('\n  # File: '+cf)

        with np.load(cf) as data:
            rdate = int( data['time'] )
            nPnts =      data['Npoints']
            if kf==0:
                corigin = str(data['origin'])
                reskm   = int(data['reskm_nmnl'])
        fb = path.basename(cf)
        vf = split('_',fb)


        cdateh = split('-',vf[-2])[0]
        list_date.append(split('h',cdateh)[0])

        kBatchName[kf] = vf[-4]

        kiDate[kf] = rdate
        kNbPoints[kf] = nPnts

        print('   * Batch: '+kBatchName[kf] )
        print('   * Date = ',mjt.epoch2clock(kiDate[kf]))
        print('   * Nb. of points = ',kNbPoints[kf] )

        kf = kf+1

    print('\n')

    nP = np.sum(kNbPoints)
    print('  ==> Total number of points:', nP)
    print('  ==> list of dates:', list_date[:])

    cdt1, cdt2 = list_date[0],list_date[-1]
    cperiod = cdt1+'-'+cdt2

    dtbin=0
    cinfobin = ''
    if cdtbin!='nemoTsi3_NoBin':
        dtbin=int(cdtbin)
        cinfobin = '_dt'+cdtbin

    if str(reskm) != creskm:
        print('ERROR: spatial scale (km) passed as argument does not match that found in deformation files!')
        exit(0)

    # Now that we know the total number of points we can allocate and fill arrays for divergence and shear
    Zshr = np.zeros(nP)
    ZDiv = np.zeros(nP)
    Ztot = np.zeros(nP)
    ZAqd = np.zeros(nP)
    Zdat = np.zeros(nP,dtype=int)

    jP = 0 ; # Counter from 0 to nP-1
    kf = 0 ; # Counter for files, 0 to nbF-1
    for cf in listnpz:
        jPe = jP + kNbPoints[kf]
        with np.load(cf) as data:
            zdiv  =      data['divergence']
            zshr  =      data['shear']
            ztot  =      data['total']
            za    =      data['quadArea']
        #kf
        ZDiv[jP:jPe] = cfg.rc_day2sec*zdiv ; # day^-1
        Zshr[jP:jPe] = cfg.rc_day2sec*zshr ; # day^-1
        Ztot[jP:jPe] = cfg.rc_day2sec*ztot ; # day^-1
        ZAqd[jP:jPe] =     za     ; # km^2
        Zdat[jP:jPe] =  kiDate[kf]
        #
        jP = jPe
        kf = kf+1

    del zdiv, zshr, ztot, za


    cfroot = corigin+cinfobin+'_'+str(reskm)+'km_'+cperiod


    # Save it for the scaling and PDFs to come ...
    print(' *** Saving: '+cd_in+'/def_SHR_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_SHR_'+cfroot+'.npz', name='shear', origin=corigin,
                         Nbatch=nbF, dates_batch=kiDate,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nP, xshr=Zshr, quadArea=ZAqd, dates_point=Zdat )

    print(' *** Saving: '+cd_in+'/def_DIV_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_DIV_'+cfroot+'.npz', name='divergence', origin=corigin,
                         Nbatch=nbF, dates_batch=kiDate,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nP, xdiv=ZDiv, quadArea=ZAqd, dates_point=Zdat )

    print(' *** Saving: '+cd_in+'/def_TOT_'+cfroot+'.npz')
    np.savez_compressed( cd_in+'/def_TOT_'+cfroot+'.npz', name='shear', origin=corigin,
                         Nbatch=nbF, dates_batch=kiDate,
                         reskm_nmnl=reskm, period=cperiod, dtbin=dtbin,
                         Np=nP, xtot=Ztot, quadArea=ZAqd, dates_point=Zdat )

