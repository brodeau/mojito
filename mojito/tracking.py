
import numpy as np

rmin_conc = 0.1 ; # ice concentration below which we disregard the point...

rFoundKM = 2.5

def find_ji_of_min(x):
    '''
    # Yes, reinventing the wheel here, but it turns out
    # it is faster this way!
    '''
    k = x.argmin()
    nx = x.shape[1]
    return k//nx, k%nx



def GetTimeSpan( dt, vtime_mod, iSdA, iMdA, iMdB, iStop=None, iverbose=0 ):
    '''
         * iStop: date at which to stop !
                  => if not provided we stop at latest model date!!!
    '''
    ltStop = ( iStop )
    #    
    # Confrontation of time info seeding and model:
    from .util import epoch2clock
    #
    if iSdA < iMdA-dt/2 or iSdA> iMdB-dt/2:
        print('PROBLEM: time in the seeding file ('+epoch2clock(iSdA)+') is outside of what model spans!')
        exit(0)
    #
    kt0 = np.argmin(np.abs(vtime_mod[:]-iSdA)) + 1
    print('    * [GetTimeSpan]: will start using model record',kt0,'of SI3 file =>',epoch2clock(vtime_mod[kt0]))
    #
    if ltStop:
        ktN = np.argmin(np.abs(vtime_mod[:]-iStop))
    else:
        ktN = len(vtime_mod) - 1        
    print('    * [GetTimeSpan]: will stop at record',ktN,' =>',epoch2clock(vtime_mod[ktN]))
    Nt = ktN - kt0 + 1
    print('       ==> '+str(Nt)+' model records')
    print('       ==> that makes '+str((vtime_mod[ktN]-vtime_mod[kt0])/(3600.*24))+' days of ice particule tracking.')
    #
    return Nt, kt0, ktN



def IsInsideQuadrangle( y, x, quad ):
    '''
        * quad: shape(4,2) [ [y0,x0], [y1,x1], [y2,x2], [y3,x3] ]    
    '''
    n = len(quad)
    if len(quad) != 4:
        print('ERROR: `len(quad) !=: 4`') ; exit(0)
    #
    lInside = False
    p2x = 0.0
    z2y = 0.0
    xints = 0.0
    [z1y,z1x] = quad[0,:]
    #
    for i in range(n+1):
        z2y,z2x = quad[i%n,:]
        if y > min(z1y,z2y):
            if y <= max(z1y,z2y):
                if x <= max(z1x,z2x):
                    if z1y != z2y:
                        xints = (y-z1y)*(z2x-z1x)/(z2y-z1y) + z1x
                    if z1x == z2x or x <= xints:
                        lInside = not lInside
        #
        z1y, z1x = z2y, z2x
        #
    return lInside





# Future fancy "Find containing cell"!


def TheCell( pyx, kjiT, pYf, pXf, iverbose=0 ):
    ''' 
        # LOLO: I want it to use Lat,Lon rather than Y,X !!!

        Based on the nearest T-point, find the cell/mesh (defined as the polygon that
        joins 4 neighbor F-points) that contains the target point.
        
        Mind: the indices of the T-point we return is that of the center of the
              identified cell! And in very unusual cases, this is not the same
              as the nearest T-point provided as an input...
       Input:
        * pyx      : (y,x) target point cartesian coordinates [km]
        * kjiT     : (j,i) indices of nearest T-point to (y,x) that was found...
        * pYf, pXf : 2D arrays of cartesian coordinates of the model F-point grid
       Output:
        * lPin    : success or not (boolean)
        * [jT,iT] : actual T-point at the center of the identified cell
        * [[jf1,jf2,jf3,jf4],[if1,if2,if3,if4]]: the 4 F-points that make the vertices
                                                 of the identified cell (ConterClockwize starting from BLC)
    '''
    (zy,zx) = pyx
    (kj,ki) = kjiT
    #
    lPin = False
    kp=0             ; # pass counter
    while (not lPin) and (kp<5):
        kp = kp + 1
        if iverbose>0 and kp>1: print('  * [SeedInit()] search for proper F-cell => test option #'+str(kp)+'!')
        if   kp==1:
            jT,iT =kj,ki   ; # Normal case!!! 99% of all cases !!!
        elif kp==2:
            jT,iT =kj,ki+1 ; # maybe F-cell is the one next to the right?
        elif kp==3:
            jT,iT =kj+1,ki ; # maybe F-cell is the one next above?
        elif kp==4:
            jT,iT =kj,ki-1 ; # maybe F-cell is the one next to the left?
        elif kp==5:
            jT,iT =kj-1,ki ; # maybe F-cell is the one next below?
        #
        # Based on this @center T-point (here `jT,iT`), the proper cell/mesh
        # should be the polygon defined by the 4 following F-points:
        # (indexing is anti-clockwize, starting from bottom left F-point)
        [jf1,jf2,jf3,jf4] = [ jT-1, jT-1, jT, jT  ]
        [if1,if2,if3,if4] = [ iT-1, iT,   iT, iT-1 ]                    
        #
        zquad = np.array( [ [pYf[jf1,if1],pXf[jf1,if1]], [pYf[jf2,if2],pXf[jf2,if2]],
                            [pYf[jf3,if3],pXf[jf3,if3]], [pYf[jf4,if4],pXf[jf4,if4]] ] )
        lPin = IsInsideQuadrangle( zy, zx, zquad )
        #
    if iverbose>0:
        if kp>1 and lPin: print('        => option #'+str(kp)+' did work!  :)')
    #
    return lPin, [jT,iT], [[jf1,jf2,jf3,jf4],[if1,if2,if3,if4]]





# This is to become a generic "Find containing Cell" for NEMO...
# => can be used to find the T-centric cell (with 4 vertices being F-points)
# => or to find the F-centric cell (with 4 vertices being T-points) => what we need for linear interpolation of a T-field onto point inside this cell...
# CALL:      pntID=pIDs[jP]
def FCC( pntGcoor, pLat, pLon, pLatC, pLonC, cellType='T', rd_found_km=10., resolkm=[],
         ji_prv=(), np_box_r=10, max_itr=5, pntID=None, iverbose=0 ):
    '''
        Provided an input target-point geographic coordinates, locate the grid mesh cell that
        includes this point.
        The grid mesh is of type Arakawa C-grid, typically a NEMO/ORCA type of grid...

    INPUT:
             * pntGcoor  : GPS coordinates (lat,lon) of target point    ([real],[real])
             * pLat      : array of source grid lat.  @T/F-points if `cellType='T/F'`) [real]
             * pLon      : array of source grid long. @T/F-points if `cellType='T/F'`) [real]
             * pLatC     : array of source grid lat.  @F/T-points if `cellType='T/F'`) [real]
             * pLonC     : array of source grid long. @F/T-points if `cellType='T/F'`) [real]
             * resolkm   : array of source grid approximate local resolution [km] [real]
                           because grids like ORCA of NEMO can have strong spatial heterogenity of resolution...
    RETURNS:
             * j,i : indices of the grid mesh cell containing the target point
                     => -1,-1 if something went wrong or point was not found
    '''
    (zlat,zlon) = pntGcoor
    
    lbla = ( iverbose>0 and pntID )
    
    cP0 = cellType      ; # string for type of center point 
    if cellType=='T':
        cP4C = 'F'      ; # string for type of corner/vertices points
    elif cellType=='F':
        cP4C = 'T'      ; # string for type of corner/vertcices points
    else:
        print('ERROR [FCC]: for now we just expect the mesh to be centered on "T" or "F" points.')
        exit(0)

    icancel = 0
    
    # First, find the nearest `cellType`-point to `pntGcoor`:
    (jX, iX) = NearestPoint( pntGcoor, pLat, pLon, rd_found_km=rd_found_km, resolkm=resolkm,
                             ji_prv=ji_prv, np_box_r=np_box_r, max_itr=max_itr )

    if jX>=0 and iX>=0:
        # Ok a nearest point was found!                
        if iverbose>0:
            print('     ==> nearest '+cP0+'-point for ',zlat,zlon,' on target grid:', jX, iX, '==> lat,lon:',
                  round(pLat[jX,iX],3), round(pLon[jX,iX],3))

        # Here problem is to have `TheCell` working using Lat,Lon rather than projected Y,X as it is done in `FindContainingCell` !!!
        # => we apply a stereographic projection of a litle region centered on the nearest point!

        
        exit(0)
        lPin, [jM,iM], [[jc1,jc2,jc3,jc4],[ic1,ic2,ic3,ic4]] = TheCell( (zlat,zlon), (jX,iX), pLatC, pLonC, iverbose=iverbose )
        #
        if not lPin:
            print('WARNING [SeedInit()]: could not find the proper F-point cell!!!')
            print('         => when lookin for point:',zlat,zlon)
            kcancel[jP] = 0
            if iverbose>0: print('        ===> I CANCEL buoy '+str(pntID)+'!!! (NO proper F-point cell found)')


            


    else:
        if lbla: print('  [FCC]      ===> I CANCEL buoy '+str(pntID)+'!!! (NO nearest '+cP0+'-point found for ',zlat,zlon,')')




        
                        
    #if kcancel[jP] == 1:            
    #    # Tests for canceling buoy or not:
    #    icncl = Survive( pntID, [jX,iX] , maskT, pIceC=xIceConc, iverbose=iverbose )
    #    if icncl>0: kcancel[jP] = 0
        
    #if kcancel[jP] == 1:
        # Everything is okay, now we locate the cell/mesh (polygon joining 4 F-points) that includes
        # our target point (zy,zx)





    #
    return 0
    



def NearestPoint( pntGcoor, pLat, pLon, rd_found_km=10., resolkm=[], ji_prv=(), np_box_r=10, max_itr=5 ):
    '''
    # * pntGcoor : GPS coordinates (lat,lon) of target point    ([real],[real])
    # * pLat        : array of source grid latitude            2D numpy.array [real]
    # * pLon        : array of source grid longitude           2D numpy.array [real]
    # * resolkm   : array of source grid approximate local resolution [km] 2D numpy.array [real]
    #               because grids like ORCA of NEMO can have strong spatial heterogenity of resolution...
    '''
    from .util import Haversine
    #
    (Ny,Nx) = pLat.shape
    if np.shape(pLon) != (Ny,Nx):
        print('ERROR [NearestPoint]: `pLat` & `pLon` do not have the same shape!')
        exit(0)
    l2Dresol = ( np.shape(resolkm)==(Ny,Nx) )
    lbox     = ( len(ji_prv)==2 )
    #
    (latP,lonP) = pntGcoor
    #        
    if lbox:
        (j_prv,i_prv) = ji_prv
        j1, j2 = max(j_prv-np_box_r,0), min(j_prv+np_box_r+1,Ny)
        i1, i2 = max(i_prv-np_box_r,0), min(i_prv+np_box_r+1,Nx)
    else:
        (j1,i1 , j2,i2) = (0,0 , Ny,Nx)
    #
    jy, jx = -1,-1 ; # "not found" flag value...
    lfound = False    
    rfnd   = rd_found_km
    igo    = 0
    #
    while (not lfound) and igo<max_itr :
        igo = igo + 1
        if lbox and igo>1:
            (j1,i1 , j2,i2) = (0,0 , Ny,Nx) ; # Falling back on whole domain for second pass...
        xd = Haversine( latP, lonP,  pLat[j1:j2,i1:i2], pLon[j1:j2,i1:i2] )
        jy, jx = find_ji_of_min( xd )
        #
        if igo==1 and l2Dresol: rfnd = 0.5*resolkm[jy,jx]
        #
        if not lbox and igo==1: igo=2 ; # we jump one round because no need of the pass to global...
        #
        lfound = ( xd[jy,jx] < rfnd )
        if igo>1 and not lfound:
            rfnd = 1.2*rfnd ; # increasing validation distance criterion by 20 %
    #
    if lbox:
        jy, jx = jy+j1, jx+i1 ; # found in the zoom box => translate to indices in whole domain:
    #
    if jy<0 or jx<0 or jy>=Ny or jx>=Nx or igo==max_itr:
        if ivrb>0: print('    WARNING [NearestPoint()]: did not find a nearest point for target point ',latP,lonP,' !')
        if ivrb>0: print('            => last tested distance criterions =', rfnd,' km')
        jy, jx = -1,-1
    #
    return (jy, jx)


def _ccw_( pcA, pcB, pcC ):
    '''
        * pcA: coordinates of point A => [y_A,x_A]
        * etc...
    '''
    return (pcC[0]-pcA[0])*(pcB[1]-pcA[1]) > (pcB[0]-pcA[0])*(pcC[1]-pcA[1])

def intersect2Seg( pcA, pcB, pcC, pcD ):
    '''
     Return true if line segments AB and CD intersect

        * pcA: coordinates of point A => [y_A,x_A]
        * etc...
    '''
    return ( _ccw_(pcA,pcC,pcD) != _ccw_(pcB,pcC,pcD) ) and ( _ccw_(pcA,pcB,pcC) != _ccw_(pcA,pcB,pcD) )



def Survive( kID, kjiT, pmskT, pIceC=[],  iverbose=0 ):
    '''
        Suite of tests to decide whether a buoy should be killed or not...
    '''
    (Nj,Ni) = np.shape(pmskT) ; # shape of the NEMO domain...
    
    [jT,iT] = kjiT
    ikill   = 0

    # Test: too close to NEMO domain boundaries:
    if ikill==0:
        if jT in [0,1,Nj-2,Nj-1] or iT in [0,1,Ni-2,Ni-1]:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (reaching NEMO domain boundaries)')
    
    # Test on land-sea mask / continent:
    if ikill==0:
        zmt = pmskT[jT,iT] + pmskT[jT,iT+1]+pmskT[jT+1,iT]+pmskT[jT,iT-1]+pmskT[jT-1,iT-1]
        if zmt < 5:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (too close to or over the land-sea mask)')

    # Test on sea-ice concentration:
    if ikill==0:
        if len(np.shape(pIceC))==2:
            zic = 0.2*(pIceC[jT,iT] + pIceC[jT,iT+1]+pIceC[jT+1,iT]+pIceC[jT,iT-1]+pIceC[jT-1,iT-1])
            if iverbose>1: print('     =>> 5P sea-ice concentration =',zic)
        if zic < rmin_conc:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (mean 5P ice concentration at T-point of the cell =',zic,')')
    #
    return ikill




def FindContainingCell( pyx, kjiT, pYf, pXf, iverbose=0 ):
    ''' 
        Based on the nearest T-point, find the cell/mesh (defined as the polygon that
        joins 4 neighbor F-points) that contains the target point.
        
        Mind: the indices of the T-point we return is that of the center of the
              identified cell! And in very unusual cases, this is not the same
              as the nearest T-point provided as an input...
       Input:
        * pyx      : (y,x) target point cartesian coordinates [km]
        * kjiT     : (j,i) indices of nearest T-point to (y,x) that was found...
        * pYf, pXf : 2D arrays of cartesian coordinates of the model F-point grid
       Output:
        * lPin    : success or not (boolean)
        * [jT,iT] : actual T-point at the center of the identified cell
        * [[jf1,jf2,jf3,jf4],[if1,if2,if3,if4]]: the 4 F-points that make the vertices
                                                 of the identified cell (ConterClockwize starting from BLC)
    '''
    (zy,zx) = pyx
    (kj,ki) = kjiT
    #
    lPin = False
    kp=0             ; # pass counter
    while (not lPin) and (kp<5):
        kp = kp + 1
        if iverbose>0 and kp>1: print('  * [SeedInit()] search for proper F-cell => test option #'+str(kp)+'!')
        if   kp==1:
            jT,iT =kj,ki   ; # Normal case!!! 99% of all cases !!!
        elif kp==2:
            jT,iT =kj,ki+1 ; # maybe F-cell is the one next to the right?
        elif kp==3:
            jT,iT =kj+1,ki ; # maybe F-cell is the one next above?
        elif kp==4:
            jT,iT =kj,ki-1 ; # maybe F-cell is the one next to the left?
        elif kp==5:
            jT,iT =kj-1,ki ; # maybe F-cell is the one next below?
        #
        # Based on this @center T-point (here `jT,iT`), the proper cell/mesh
        # should be the polygon defined by the 4 following F-points:
        # (indexing is anti-clockwize, starting from bottom left F-point)
        [jf1,jf2,jf3,jf4] = [ jT-1, jT-1, jT, jT  ]
        [if1,if2,if3,if4] = [ iT-1, iT,   iT, iT-1 ]                    
        #
        zquad = np.array( [ [pYf[jf1,if1],pXf[jf1,if1]], [pYf[jf2,if2],pXf[jf2,if2]],
                            [pYf[jf3,if3],pXf[jf3,if3]], [pYf[jf4,if4],pXf[jf4,if4]] ] )
        lPin = IsInsideQuadrangle( zy, zx, zquad )
        #
    if iverbose>0:
        if kp>1 and lPin: print('        => option #'+str(kp)+' did work!  :)')
    #
    return lPin, [jT,iT], [[jf1,jf2,jf3,jf4],[if1,if2,if3,if4]]



def SeedInit( pIDs, pSG, pSC, platT, plonT, pYf, pXf, pResolKM, maskT, xIceConc=[], iverbose=0 ):
    '''
    '''
    from time import time
    #
    (nP,n2) = np.shape(pSG)
    #
    if np.shape(pSC)!=(nP,n2):
        print('ERROR [SeedInit]: shape disagreement for `pSG` and `pSC`!'); exit(0)
    if n2!=2:
        print('ERROR [SeedInit]: wrong shape for `pSG` and `pSC`!'); exit(0)
    #
    #                  jF,iF:       # j,i indices of the F-point that defines the current mesh/cell
    #                                 # (was the nearest point when we searched for nearest point,
    #                                 #  as buoys move moves within the cell, it might not be the nearest point)
    zjiT    = np.zeros((nP,2),   dtype=int)
    zJIvrt  = np.zeros((nP,2,4), dtype=int)
    kcancel = np.zeros( nP ,     dtype='i1') + 1

    time_start = time()
    
    for jP in range(nP):

        if jP%100==0:
            print('   * [SeedInit()]: finding nearest T-point and host cell... buoy #'+str(jP)+' / '+str(nP)+' ...')
        
        if iverbose>0:
            print('   * [SeedInit()]: focus on buoy with ID:'+str(pIDs[jP]))

        # Initial position of the buoy:
        zlat,zlon  = pSG[jP,0], pSG[jP,1] ; # degrees!
        zy  ,zx    = pSC[jP,0], pSC[jP,1] ; # km !
        zic = 1.
        
        # 1/ Nearest T-point on NEMO grid:
        (jT, iT) = NearestPoint( (zlat,zlon), platT, plonT, rd_found_km=rFoundKM, resolkm=pResolKM, max_itr=10 )

        if jT<0 or iT<0:
            kcancel[jP] = 0
            if iverbose>0: print('        ===> I CANCEL buoy '+str(pIDs[jP])+'!!! (NO nearest T-point found for ',zlat,zlon,')')
                            
        if kcancel[jP] == 1:            
            # Ok a nearest point was found!    
            if iverbose>0:
                print('     ==> nearest T-point for ',zlat,zlon,' on NEMO grid:', jT, iT, '==> lat,lon:',
                      round(platT[jT,iT],3), round(plonT[jT,iT],3))
            # Tests for canceling buoy or not:
            icncl = Survive( pIDs[jP], [jT,iT] , maskT, pIceC=xIceConc, iverbose=iverbose )
            if icncl>0: kcancel[jP] = 0
            
        if kcancel[jP] == 1:
            # Everything is okay, now we locate the cell/mesh (polygon joining 4 F-points) that includes
            # our target point (zy,zx)
            lPin, zjiT[jP,:], zJIvrt[jP,:,:] = FindContainingCell( (zy,zx), (jT,iT), pYf, pXf, iverbose=iverbose )
            #
            if not lPin:
                print('WARNING [SeedInit()]: could not find the proper F-point cell!!!')
                print('         => when lookin for point:',zlat,zlon)
                kcancel[jP] = 0
                if iverbose>0: print('        ===> I CANCEL buoy '+str(pIDs[jP])+'!!! (NO proper F-point cell found)')
                        
    ### for jP in range(nP)
    time_stop = time()
    print(' * [SeedInit]: number of seconds it took to locate all the points in a target grid cell:', time_stop-time_start )
    
    # Now we can shrink the arrays based on `kcancel`
    iKeep = np.arange(nP,dtype=int)
    nPn   = np.sum(kcancel)
    if nPn<nP:
        (iGone,) = np.where(kcancel==0)
        print(' * [SeedInit()]: '+str(nP-nPn)+' "to-be-seeded" buoys have to be canceled.')
        print('          => their IDs:',pIDs[iGone])
        print('        ==> need to shrink some arrays, number of valid buoys is now',nPn)
        nP = nPn
        (iKeep,) = np.where(kcancel==1)
                
    return nP, pSG[iKeep,:], pSC[iKeep,:], pIDs[iKeep], zjiT[iKeep,:], zJIvrt[iKeep,:,:]


def CrossedEdge( pP1, pP2, ji4vert, pY, pX,  iverbose=0 ):
    '''
        * pP1    : point 1 as [y1,x1]
        * pP2    : point 2 as [y2,x2]
        * ji4vert : array(2,4), containg the 4 j,i indices of the 4 vertices of the mesh
    '''
    #
    for kk in range(4):
        kp1 = (kk+1)%4
        [j1,i1] = ji4vert[:,kk]
        [j2,i2] = ji4vert[:,kp1]
        #
        ll = intersect2Seg( pP1, pP2,  [pY[j1,i1],pX[j1,i1]], [pY[j2,i2],pX[j2,i2]] )
        if ll: break
        #
    if iverbose>0:
        vdir = ['bottom', 'right-hand', 'upper', 'left-hand']
        print('    [CrossedEdge()]: particle is crossing the '+vdir[kk]+' edge of the mesh!')
    return kk+1


def NewHostCell( kcross, pP1, pP2, ji4vert, pY, pX,  iverbose=0 ):
    '''
       Find in which adjacent cell, the point has moved into !

       * kcross : integer describing which (old) cell edge the point has crossed

       In rare cases, the buoy could possibly move into diagonally adjacent cells
       (not only bottom,right,upper,left cells...)
       These "diagonally-adjacent" meshes yields output: 5,6,7,8
    '''
    #
    [ [ jbl, jbr, jur, jul ], [ ibl, ibr, iur, iul ] ] = ji4vert[:,:]
    #
    knhc = kcross
    if  kcross==1:
        # Crosses the bottom edge:
        if   intersect2Seg( pP1, pP2, [pY[jbl,ibl],pX[jbl,ibl]], [pY[jbl-1,ibl],pX[jbl-1,ibl]] ):
            knhc=5 ; # bottom left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jbr,ibr],pX[jbr,ibr]], [pY[jbr-1,ibr],pX[jbr-1,ibr]] ):
            knhc=6 ; # bottom right diagonal
        #
    elif kcross==2:
        # Crosses the RHS edge:
        if   intersect2Seg( pP1, pP2, [pY[jbr,ibr],pX[jbr,ibr]], [pY[jbr,ibr+1],pX[jbr,ibr+1]] ):
            knhc=6 ; # bottom right diagonal
        elif intersect2Seg( pP1, pP2, [pY[jur,iur],pX[jur,iur]], [pY[jur,iur+1],pX[jur,iur+1]] ):
            knhc=7 ; # bottom right diagonal
        #
    elif kcross==3:
        # Crosses the upper edge:
        if   intersect2Seg( pP1, pP2, [pY[jul,iul],pX[jul,iul]], [pY[jul+1,iul],pX[jul+1,iul]] ):
            knhc=8 ; # upper left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jur,iur],pX[jur,iur]], [pY[jur+1,iur],pX[jur+1,iur]] ):
            knhc=7 ; # bottom right diagonal
        #
    elif kcross==4:
        # Crosses the LHS edge:
        if   intersect2Seg( pP1, pP2, [pY[jul,iul],pX[jul,iul]], [pY[jul,iul-1],pX[jul,iul-1]] ):
            knhc=8 ; # upper left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jbl,ibl],pX[jbl,ibl]], [pY[jbl,ibl-1],pX[jbl,ibl-1]] ):
            knhc=5 ; # bottom left diagonal
    #
    if iverbose>0:
        vdir = ['bottom', 'RHS', 'upper', 'LHS', 'bottom-LHS', 'bottom-RHS', 'upper-RHS', 'upper-LHS' ]
        print('    *** Particle is moving into the '+vdir[knhc-1]+' mesh !')
    #
    return knhc



def UpdtInd4NewCell( knhc, ji4vert, kjiT ):
    '''
        Update the mesh indices according to the new host cell
    '''
    if   knhc==1:
        # went down:
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[0] = kjiT[0]-1
    elif knhc==5:
        print(' * [UpdtInd4NewCell()]: WE HAVE A 5 !!!!')
        # went left + down
        ji4vert[1,:] = ji4vert[1,:] - 1
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[1] = kjiT[1]-1
        kjiT[0] = kjiT[0]-1
    elif knhc==6:
        print(' * [UpdtInd4NewCell()]: WE HAVE A 6 !!!!')
        # went right + down
        ji4vert[1,:] = ji4vert[1,:] + 1
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[1] = kjiT[1]+1
        kjiT[0] = kjiT[0]-1
    elif knhc==2:
        # went to the right
        ji4vert[1,:] = ji4vert[1,:] + 1
        kjiT[1] = kjiT[1]+1
    elif knhc==3:
        # went up:
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[0] = kjiT[0]+1
    elif knhc==7:
        print(' * [UpdtInd4NewCell()]: WE HAVE A 7 !!!!')
        # went right + up
        ji4vert[1,:] = ji4vert[1,:] + 1
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[1] = kjiT[1]+1
        kjiT[0] = kjiT[0]+1
    elif knhc==8:
        print(' * [UpdtInd4NewCell()]: WE HAVE A 8 !!!!')
        # went right + up
        ji4vert[1,:] = ji4vert[1,:] - 1
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[1] = kjiT[1]-1
        kjiT[0] = kjiT[0]+1
    elif knhc==4:
        # went to the left
        ji4vert[1,:] = ji4vert[1,:] - 1
        kjiT[1] = kjiT[1]-1
    else:
        print('ERROR: unknown direction, knhc=',knhc)
        exit(0)

    return ji4vert, kjiT





# Misc. seeding functions (when seeding not based on input file)...

def debugSeeding():
    zLatLon = np.array([
        [84. ,  20.],
        [89. ,  50.],
        [63. , -11.], # point outside of domain
        [85. , 100.],        
        [89. , 100.],
        [79. , 180.],
        [76. ,  46.], # No ice, Barents Sea
        [75. , 190.],
        [85.2, -15.], # Too close to land
        [75. , 210.],
        [75. , -72.], # Baffin bay
        [83. , 200.],
        [79. , -42.], # over mask (Greenland!)
        [85. , 300.]  ])
    return zLatLon

def debugSeeding1():
    zLatLon = np.array([ [75.,190.] ])
    return zLatLon


def nemoSeed( pmskT, platT, plonT, pIC, khss=1, fmsk_rstrct=None ):

    zmsk = pmskT[::khss,::khss]
    zlat = platT[::khss,::khss]
    zlon = plonT[::khss,::khss]

    (Nj,Ni) = np.shape(zmsk)
    ztmp = np.zeros((Nj,Ni))
    
    #if fmsk_rstrct:
    #    with Dataset(fmsk_rstrct) as id_mr:
    #        maskR = id_mr.variables['mask'][::khss,::khss]
    #        if np.shape(maskR) != (Nj,Ni):
    #            print('ERROR [nemoSeed()]: restricted area mask does not agree in shape with model output!'); exit(0)
    
    msk_T = np.zeros((Nj,Ni), dtype='i1')

    msk_T[:,:] = zmsk[:,:]

    # Only north of ...
    msk_T[np.where(zlat < 55.)] = 0
        
    # Only over a decent concentration of ice:
    ztmp[:,:] = pIC[::khss,::khss]
    msk_T[np.where(ztmp < 0.9)] = 0

    (idy_keep, idx_keep) = np.where( msk_T==1 )
    NbIPt = len(idy_keep)
    zLatLon = np.zeros((NbIPt,2))
    
    for jp in range(NbIPt):
        jj = idy_keep[jp]
        ji = idx_keep[jp]
        zLatLon[jp,:] = [ zlat[jj,ji],zlon[jj,ji] ]
        
    del zmsk,zlat,zlon
        
    return zLatLon





