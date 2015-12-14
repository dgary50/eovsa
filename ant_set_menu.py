import requests
import time

def ant_set_menu(ant,item,value):

    ant_ip = ant+'.solar.pvt'
    # Fill in your details here to be posted to the login form.
    if ant == 'ant12':
    loginfo = {
        'UserLogin_ForceIP':'0',
        'UserLogin_InputUsername': 'root',
        'UserLogin_InputPassword': 'intertronic'
    }

    # Use 'with' to ensure the session context is closed after use.
    with requests.Session() as s:
    r = s.post('http://ant12.solar.pvt/Forms/login_1', data=loginfo)
    # print the html returned or something more intelligent to see if it's a successful login page.
    lines = r.text.split('\n')
    print 'Login status:'
    for line in lines:
        if line.find('LOG') != -1:
           print line+'\n'
    loginval = r.url.split('/')[-2]
    print loginval

#    # Set menu 10.00 as 1070
#    payload = {'NewValue':'1070'}
#    url = 'http://ant12.solar.pvt/US/1000/Forms/parameter_3'
#    r = s.post(url,data=payload)
#    print 'Post Newvalue to 10.00 status:',r.status_code

#    time.sleep(2)

#    # Set menu 10.38 as 100
#    payload = {'NewValue':'100'}
#    url = 'http://ant12.solar.pvt/US/1038/Forms/parameter_3'
#    r = s.post(url,data=payload)
#    print 'Post Newvalue to 10.38 status:',r.status_code

    # Set menu 20.40 as 777
    payload = {'NewValue':'777'}
    url = 'http://ant12.solar.pvt/US/'+loginval+'/2040/Forms/parameter_3'
    print url
    r = s.post(url,data=payload)
    lines = r.text.split('\n')
    print 'Post Newvalue to 20.40:'
    for line in lines:
        if line.find('<br>') != -1:
            print line+'\n'

#http://ant12.solar.pvt/Forms/login_1?UserLogin_ForceIP=0&UserLogin_InputUsername=root&UserLogin_InputPassword=intertronic
#http://ant12.solar.pvt/US/2039/Forms/parameter_3?NewValue=712


