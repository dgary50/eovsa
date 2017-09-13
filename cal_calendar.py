def find_calibrations(year, month):
    import calendar
    from util import Time
    import cal_header as ch
    import stateframe
    hc = calendar.HTMLCalendar(calendar.SUNDAY)
    html_table = hc.formatmonth(year, month)
    lines = html_table.split('\n')
    lines[2] = lines[2].replace('th class','th width="100" class')
    c = calendar.TextCalendar(calendar.SUNDAY)
    for i in c.itermonthdays(year,month):
        if i != 0:
            t = Time(str(year)+'-'+str(month)+'-'+str(i)+' 20:00')
            cals = []
            for caltype in [8,9,10]:
                xml, buf = ch.read_cal(caltype,t)
                if buf is None:
                    cals.append(0)
                else:
                    tout = Time(stateframe.extract(buf,xml['Timestamp']),format='lv')
                    if (t - tout).value < 1./3:
                        cals.append(1)
                    else:
                        cals.append(0)
            for k,line in enumerate(lines[3:]):
                idx = line.find(str(i))
                ns = len(str(i))
                if idx != -1:
                    line = line[:idx+ns]+'<br>-r- -p- -tp- <br>&nbsp;{} &nbsp; {} &nbsp; {}'.format(*cals)+line[idx+ns:]
                    break
            lines[k+3] = line
    #print ''.join(line+'\n' for line in lines)
    return ''.join(line+'\n' for line in lines)

if __name__ == '__main__':
    import sys
    import numpy as np
    try:
        year = np.int(sys.argv[1])
        month = np.int(sys.argv[2])
    except:
        print 'Error interpreting command line argument'
    txt = find_calibrations(year,month)
    f = open('/common/webplots/cal_status/{:4}{:02d}'.format(year,month)+'.txt','w')
    f.write(txt)
    f.close()