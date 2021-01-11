def test_file(file):
    import aipy
    #file = '/dppdata1/IDB/IDB20170515010201'
    #file = '/dppdata1/IDB/IDB20170519010019'
    uv = aipy.miriad.UV(file)
    utcount = 0
    ut = 0.0
    try: 
        for preamble, data in uv.all():
            # Look for time change
            if preamble[1] != ut:
                ut = preamble[1]
                utcount = utcount+1
                print utcount
            #endif
        #endfor
    except:
        print "oops, file didn't work"
    #end
#END

