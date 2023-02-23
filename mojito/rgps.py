

def streamSummaryRGPS( pNBini, pTcini, pIDs, pNRc ):
    from numpy import max, shape
    from climporn import epoch2clock
    #
    (Nstrm,NbMax) = shape(pIDs)
    if Nstrm != len(pNBini):
        print('ERROR: [streamSummaryRGPS()] => error #1')
    if not NbMax == max(pNBini):
        print('ERROR: [streamSummaryRGPS()] => error #2')
    print('\n ==========   SUMMARY   ==========')
    print(' *** Number of identified streams: '+str(Nstrm))
    print(' *** Number of buoys selected in each stream:')
    for js in range(Nstrm):
        cTc0  = epoch2clock(pTcini[js])
        print('        * Stream #'+str(js)+' initiated at time bin centered around '+cTc0+' => has '+str(pNBini[js])+' buoys')
    print(' *** Max number of buoys possibly found in a stream = ',NbMax)
    print('     * shape of ZIDs =', shape(pIDs))
    print('     * shape of ZNRc =', shape(pNRc))
    print(' ===================================\n')

