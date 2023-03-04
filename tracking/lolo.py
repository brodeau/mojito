
def GetTimeSpan( dt, vtime_mod, iSdA, iSdB, iMdA, iMdB, iverbose=0 ):
    # Confrontation of time info seeding and model:
    if iSdA < iMdA-dt/2 or iSdA> iMdB-dt/2:
        print('PROBLEM: time in the seeding file ('+mjt.epoch2clock(iSdA)+') is outside of what model spans!')
        exit(0)
    #
    kt0 = np.argmin(np.abs(vtime_mod[:]-iSdA)) + 1
    print('    * [GetTimeSpan]: will start using record',kt0,'of SI3 file =>',mjt.epoch2clock(vtime_mod[kt0]))
    ktN = np.argmin(np.abs(vtime_mod[:]-iSdB))
    print('    * [GetTimeSpan]: will stop at record',ktN,' =>',mjt.epoch2clock(vtime_mod[ktN]))
    Nt = ktN - kt0 + 1
    print('       ==> '+str(Nt)+' model records')
    print('       ==> that makes '+str((vtime_mod[ktN]-vtime_mod[kt0])/(3600.*24))+' days of ice particule tracking.')
    #
    return Nt, kt0, ktN

