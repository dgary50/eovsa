# Name:
#   DAILY_XSP_BATCH
# Create daily "all_day" dynamic spectrum plots of EOVSA solar 
# cross-correlation data
#
#  2019-02-19 DG
#    Batch version written
#

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import os
    try:
        user_paths = os.environ['PYTHONPATH'].split(os.pathsep)
    except:
        user_paths = []
    print user_paths

def last_date():
    import glob
    from util import Time
    files = glob.glob('/common/webplots/flaremon/daily/2017/*.png')
    files.sort()
    date = '2016-12-31 20:00:00'
    if files != []:
        datstr = files[-1].split('XSP')[1]
        date = datstr[:4]+'-'+datstr[4:6]+'-'+datstr[6:8]
    return Time(date)
    
if __name__ == "__main__":
    ''' For non-interactive use, use a backend that does not require a display
        Usage python /common/python/current/daily_xsp.py <date>
        
        If optional argument date is given (as YYYY-MM-DD), data are processed for 
           that date only.
        If the date is omitted, data are processed for the previous two UT days
          (yesterday and day-before-yesterday)
    '''
    import glob, sys
    import daily_xsp
    from util import Time
    t = last_date()
    blah = {}
    while blah == {}:
        mjd = t.mjd
        t = Time(mjd+1,format='mjd')
        print 'Working on:',t.iso[:19]
        blah = daily_xsp.allday_udb(t=t, savfig=True)   # Process time t
        
