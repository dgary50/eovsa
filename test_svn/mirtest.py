import os

def valid_miriad_dset0(filelist0):
    '''Returns True or False for valid or invalid Miriad datasets,
    checks for existnce of the directory, and then for flags, header,
    vartable, and visdata. Also returns names of valid datasets'''

    if len(filelist0) == 0:
        print 'valid_miriad_dset0: No files input'
        return False
    #endif

    #need a list input, otherwise all sorts of things are messed up
    if not isinstance(filelist0, list):
        filelist = [filelist0]
    else:
        filelist = filelist0
    #endelse

    n = len(filelist)
    otp = []
    ok_filelist = []
    bad_filelist = []
    for j in range(n):
        filename = filelist[j]
        tempvar = True
        if os.path.isdir(filename) == False or os.path.isfile(filename+'/flags') == False or os.path.isfile(filename+'/header') == False or os.path.isfile(filename+'/vartable') == False or os.path.isfile(filename+'/visdata') == False:
            tempvar = False
        #end if
        otp.append(tempvar)
        if tempvar == True:
            ok_filelist.append(filelist[j])
        #end if
        if tempvar == False:
            bad_filelist.append(filelist[j])
        #end if
    #endfor
    return otp, ok_filelist, bad_filelist
#End of valid_miriad_dset0
