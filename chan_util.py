global nschan, ifbw, nschanx, nsavg
nschan = 4096
ifbw = 600.
nschanx = 341
nsavg = [64,97,131,162,200,227,262]+[284]*27


def chan_asmt(bnd):
    # Input band (bnd) ranges from 1-34, but the code below uses band,
    # which ranges from 0-33.
    band = bnd - 1
    if bnd < 1 or bnd > 34:
        # Error in band number provided
        return -1

    global nschan, ifbw, nschanx, nsavg
    df = ifbw/nschan

    chasmt = [0]*nschanx

    # Create list of science channels for each band up to this one
    nscichan = []
    for n in nsavg[0:band+1]:
        nscichan.append(int(500/(n*df)))
    # Sum number of science channels for all bands up to but
    # not including this one.
    nsci = 0
    if band is 0:
        pass
    else:
        for n in nscichan[0:band]:
            nsci += n
    for i in range(nscichan[band]):
        chasmt += [i+nsci+1]*nsavg[band]
    nrest = 4096 - len(chasmt)
    chasmt += [0]*nrest

    return chasmt

def start_freq(band):
    # Input band (bnd) ranges from 1-34
    if band < 1 or band > 34:
        # Error in band number provided
        return -1

    global nschan, ifbw, nschanx, nsavg
    df = ifbw/nschan

    sf = []

    # Create list of frequencies for each science channel in this band
    nscichan = int(500/(nsavg[band-1]*df))
    for n in range(nscichan):
        sf.append(0.450 + band*0.5 + (nschanx + nsavg[band-1]*n)*df/1000.)

    return sf

def sci_bw(band):
    # Input band (bnd) ranges from 1-34
    if band < 1 or band > 34:
        # Error in band number provided
        return -1

    global nschan, ifbw, nschanx, nsavg
    df = ifbw/nschan

    nscichan = int(500/(nsavg[band-1]*df))
    scibw = [nsavg[band-1]*df/1000.]*nscichan

    return scibw


