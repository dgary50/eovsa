import os, subprocess

def plate_solve(imgdir,obsdate):
    os.chdir('C:/Documents and Settings/Dale/My Documents/Home/2m_star_pointing')
    # Set environment variable used by pinpt_test.vbs script
    os.environ['FITSPATH'] = "C:/Documents and Settings/Dale/My Documents/Home/LX200GPS/DGary/"+imgdir
    sdir = "C:/Documents and Settings/Dale/My Documents/Dropbox/PythonCode/Star Pointing/"
    f = open(sdir+'startable-'+obsdate+'.txt','r')
    lines = f.readlines()
    f.close()
    o = open(sdir+'starsolutions-'+obsdate+'.txt','w')
    lines = lines[2:]  # Remove two header lines
    for line in lines:
        if line is '\n':
            # In case there is a blank line at end
            break
        # Split line into parts demarcated by white space
        tokens = line[17:].split(' ')
        ntok = len(tokens)
        rah, ram, ras = tokens[0].split(':')
        if tokens[2] is '':
            decd, decm, decs = tokens[3].split(':')
        else:
            decd, decm, decs = tokens[2].split(':')
        imgnum = tokens[ntok-1].strip()
        ra = float(rah) + float(ram)/60. + float(ras)/3600.
        sgn = +1.
        decd = float(decd)
        if decd < 0:
            sgn = -1.
            decd = -decd
        dec = sgn*(decd + float(decm)/60. + float(decs)/3600.)
        name = "star-"+imgnum+".fit"
        result = subprocess.Popen(["cscript","pinpt_test.vbs",name,str(ra),str(dec)],
                                  stdout=subprocess.PIPE,shell=True)
        try:
            res, err = result.communicate()
        except:
            res = result.communicate()
        out = res.split('**')
        if len(out) is 1:
            # Solution failed
            lineout = line[:43]+' Solution failed'
        else:
            # We have a solution
            ra_out, dec_out, sdate, stime, ampm = out[1].strip().split(' ')
            if len(stime) is 7:
                # If less than 10 hours, add leading zero to stime
                stime = '0'+stime
            if ampm is 'PM':
                # Time is PM, so add 12 hours to the start hour
                sth, stm, sts = stime.split(':')
                sth = str(int(sth)+12)
                stime = sth+':'+stm+':'+sts
            ra_out = float(ra_out)
            dec_out = float(dec_out)
            rah = int(ra_out)
            ram = int((ra_out - rah)*60.)
            ras = ((ra_out - rah)*60. - ram)*60.
            decd = int(dec_out)
            decm = int((dec_out - decd)*60.)
            decs = ((dec_out - decd)*60. - decm)*60.
            if dec_out > 0:
                lineout = line[:43]+'{:02d}:{:02d}:{:06.3f}  {:02d}:{:02d}:{:06.3f}'.format(rah,ram,ras,decd,decm,decs) + '  ' + stime
            else:
                lineout = line[:43]+'{:02d}:{:02d}:{:06.3f} {:03d}:{:02d}:{:06.3f}'.format(rah,ram,ras,decd,-decm,-decs) + '  ' + stime
        print lineout
        o.write(lineout+'\n')

    o.close()
    os.chdir(sdir)
