#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, makedirs
from glob import glob
import numpy as np
from re import split
#
import mojito   as mjt
from mojito import config as cfg

idebug=1
iplot=1 ; NameArcticProj='SmallArctic'

cprefixIn='DEFORMATIONS_' ; # Prefix of deformation files...

lExportCloud = False

#FracOverlap = 0.05
#FracOverlap = 0.2 ; # for 320km ?
FracOverlap = 0.7 ; # for 320km ? #best

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
    return zP1.intersects(zP2), zP1.intersection(zP2).area



def lOverlapAnyOtherQuad( pQxy, pMQ, scale=320, frac_overlap=0.1, iverbose=0 ):
    '''
       * pQxy: cartesian coordinates of 1 quad, shape = (4,2)
       * pMQ:  cartesian coordinates of a groupd of quads, shape = (nq,4,2)

       => will return true if the quad `pQxy` overlaps partially (or completely) at least 1 quad of `pMQ` !
          ==> to be considered an overlap the intersection area of the 2 quads must be > frac_overlap*(scale*scale)

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
    znmnl_area = float(scale*scale)
    #
    lovrlp = False
    for jq in range(nq):
        li, zAi = lIntersect2Quads( pQxy, pMQ[jq,:,:] )
        if li:
            if zAi<0.:
                print('ERROR [lOverlapAnyOtherQuad()]: negative intersection area!'); exit(0)
            if zAi/znmnl_area > frac_overlap:
                lovrlp = True
                break
    #
    return lovrlp




if __name__ == '__main__':

    kk = cfg.initialize()

    if not len(argv) in [6]:
        print('Usage: '+argv[0]+' <directory_input_npz_files> <SXXX> <dtbin_h> <creskm> <string_id_origin>')
        exit(0)

    cd_in  = argv[1]
    cbatch = argv[2]
    cdtbin = argv[3]
    creskm = argv[4]
    cidorg = argv[5]

    #cbatch = 'S%3.3i'%(int(kbatch))


    # Populating files and ordering them according to SS resolution:
    
    listnpz = np.array( glob(cd_in+'/'+cprefixIn+'*'+cidorg+'_'+cbatch+'_'+cdtbin+'*_*-'+creskm+'km.npz') )
    ccr = [ split('_',path.basename(cf))[-1] for cf in listnpz ] ; # list of resolution suffixes... (for odering)
    listnpz = listnpz[np.argsort(ccr)]
    # Polpulating deformation files available:
    nbF = len(listnpz)
    print('\n *** We found '+str(nbF)+' deformation files into '+cd_in+' !')
    print('      => files to gather into a single one:',listnpz,'\n')
    
    if lExportCloud:
        # nc/PointsOfQuadsOfDEF_RGPS_S006_dt72_19970119_19970122_275-320km.nc        
        listnc = np.array( glob('./nc/PointsOfQuadsOfDEF_'+cidorg+'_'+cbatch+'_'+cdtbin+'_*_*-'+creskm+'km.nc')  )
        ccr = [ split('_',path.basename(cf))[-1] for cf in listnc ] ; # list of resolution suffixes... (for odering)
        listnc = listnc[np.argsort(ccr)]
        if len(listnc) != nbF:
            print('ERROR: `len(listnc) != nbF`')

    
    # First, populate the number of deformation quadrangles in each file:
    kNbPoints = np.zeros(nbF, dtype=int ) ; # number of points in file
    vsubRes   = np.zeros(nbF, dtype=int ) ; # coarsifying resolution


    kf = 0
    for cf in listnpz:        
        css   = split( '-', split( '_', split('\.',path.basename(cf))[0] )[-1] )[0]
        vsubRes[kf] = int(css)
        csskm = css+'km'

        with np.load(cf) as data:
            nPnts =      data['Npoints']
            if kf==0:
                corigin = str(data['origin'])
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


    if str(reskm) != creskm:
        print('ERROR: `str(reskm) != creskm` !'); exit(0)

    cf_keep = str.replace( listnpz[kw], '_'+str(vsubRes[kw])+'-'+creskm, '_SLCT'+creskm)
    
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
                if not lOverlapAnyOtherQuad( zQxy, zInUseXY[:knxt,:,:], scale=reskm, frac_overlap=FracOverlap ):
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

    print('\n *** Saving new file: '+cf_keep)
    np.savez_compressed( cf_keep, time=itimeC, date=ctimeC, Npoints=nQtot,
                         Xc=zXc, Yc=zYc, X4=zX4, Y4=zY4, divergence=zdiv, shear=zshr, total=ztot,
                         quadArea=zAq, origin=corigin, reskm_nmnl=reskm )
    print('')



    if lExportCloud:
        jf = 0
        jq = 0
        for cf in listnc[idxF2K]:
            print('    - '+cf+' ('+str(vsubRes[idxF2K[jf]])+'km)'); # => indices of quads to keep:',idxF_Q2K[jf])
            ikeep = idxF_Q2K[jf]
            nQr = len(ikeep)
            jqe = jq+nQr
            #
            idate, zIDs, zGC, zXY, ztim = mjt.LoadNCdataMJT( cf, krec=0, lGetTimePos=True, convention='F' )

            print('shape zXY =>',np.shape(zXY))

            # Coordinates associated to the Quads we keep
            z4x = zX4[jq:jqe,:]
            z4y = zY4[jq:jqe,:]

            print('shape z4x =>',np.shape(z4x))

            # FINISH ME !!!
            exit(0)

            
            #
            jf += 1
            jq = jqe
        ### cf in listnpz[idxF2K]
        print('')


        #makedirs( './nc', exist_ok=True )    
        #kk = mjt.ncSaveCloudBuoys( './nc/PointsOfQuadsOfDEF_'+cfnm2+'.nc', vtim, zPids1, zPXY[:,:,1], zPXY[:,:,0], zPGC[:,:,1], zPGC[:,:,0],
        #                           xtime=xtim, fillVal=mjt.FillValue, corigin='RGPS' )



        
        exit(0)
        



    

    # Some plots:
    if iplot>0:

        cfnm = str.replace( path.basename(cf_keep), 'DEFORMATIONS_', '')
        cfnm = str.replace( cfnm, '.npz', '')
        
        zoom = 1
        
        k2 = cfg.updateConfig4Scale( reskm, mode='rgps' )

        cresinfo = '('+str(reskm)+' km)'
        
        cfdir = './figs/deformation/'+str(reskm)+'km'
        makedirs( cfdir, exist_ok=True )

        if corigin != 'RGPS':
            corigin = str.replace( corigin,'NEMO-','')
            corigin = str.replace( corigin,'_NANUK4_','-')
            
        cresinfo = '('+str(reskm)+' km)'

        nmproj=NameArcticProj


        zrx = [ np.min(zX4)-25. , np.max(zX4)+25. ]
        zry = [ np.min(zY4)-25. , np.max(zY4)+25. ]

        
        # Filled quads projected on the Arctic map:
        mjt.ShowDefQuadGeoArctic( zX4, zY4, cfg.rc_day2sec*zdiv, cfig=cfdir+'/map_zd_'+cfnm+'_Divergence.png', nmproj=NameArcticProj, cwhat='div',
                                  pFmin=-cfg.rc_div_max_fig, pFmax=cfg.rc_div_max_fig, zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': divergence '+str(cresinfo), idate=itimeC )

        mjt.ShowDefQuadGeoArctic( zX4, zY4, cfg.rc_day2sec*zshr, cfig=cfdir+'/map_zs_'+cfnm+'_Shear.png',      nmproj=NameArcticProj, cwhat='shr',
                                  pFmin=0.,      pFmax=cfg.rc_shr_max_fig,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': shear '+cresinfo, idate=itimeC )


        cix = ''
        cix ='max%5.5i'%(int( round( 10000.*np.max(cfg.rc_day2sec*ztot) , 0) ))+'_'

        mjt.ShowDefQuadGeoArctic( zX4, zY4, cfg.rc_day2sec*ztot, cfig=cfdir+'/map_zt_'+cix+cfnm+'_Total.png', nmproj=NameArcticProj, cwhat='tot',
                                  pFmin=0.,      pFmax=cfg.rc_tot_max_fig,  zoom=zoom, rangeX=zrx, rangeY=zry, unit=r'day$^{-1}$',
                                  title=corigin+': total deformation '+cresinfo, idate=itimeC )

