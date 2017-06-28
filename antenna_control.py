def get_trip(tripcode):
    ''' Given a numerical antenna controller Trip code, looks up the meaning of the
        Trip and returns documentation (description and options for solving) based
        on the Control Techniques documentation. '''
    import os
    # Open CT_Trip_Info text file (should be in PYTHONPATH)
    for p in os.environ['PYTHONPATH'].split(os.pathsep):
        if p == '': p = '.'
        try:
            f = open(p+os.sep+'CT_Trip_Info.txt','r')
            lines = f.readlines()
            f.close()
        except:
            pass
    bl = []
    for i,line in enumerate(lines):
        if line.strip() == '':
            bl.append(i)

    found = False
    for i,bline in enumerate(bl[:-3]):
        codelist = lines[bline+2].strip().split(',')
        for j in range(len(codelist)):
            if codelist[j] == str(tripcode):
                found=True
                break
        if found: break

    clues = []
    if found:
        # bline is index of the blank line above the found tripcode, 
        # so bline+1 is the start of the trip info
        trip = lines[bline+1].strip()
        desc = lines[bline+3].strip()
        n = 4
        while lines[bline+n][0] == '*':
            clues.append(lines[bline+n].strip())
            n += 1
        return {'Trip':trip,'Tripcode':tripcode,'Description':desc,'Clues':clues}
    else:
        return {'Trip':'No such tripcode','Tripcode':tripcode,'Description':'Tripcode not found','Clues':[]}
