import socket, sys, struct, time
host  = 'geobrickanta.solar.pvt'
port = 1025

rq_type = {'upload': '\xc0', 'download': '\x40'}
rq = {'sendline':    '\xb0',
      'getline':     '\xb1',
      'flush':       '\xb3',
      'getmem':      '\xb4',
      'setmem':      '\xb5',
      'setbit':      '\xba',
      'setbits':     '\xbb',
      'port':        '\xbe',
      'getresponse': '\xbf',
      'readready':   '\xc2',
      'response':    '\xc4',
      'getbuffer':   '\xc5',
      'writebuffer': '\xc6',
      'writeerror':  '\xc7',
      'fwdownload':  '\xcb',
      'ipaddress':   '\xe0'}

def send_str(msg):

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) # TCP
    s.connect((host,port))
    # Create the ethernet packet contents
    buf = rq_type['download']+rq['sendline']
    buf += struct.pack('H',0)  # Value = 0
    buf += struct.pack('H',0)  # Index = 0
    buf += struct.pack('H',socket.htons(len(msg)+1))  # Length including null termination
    buf += struct.pack(str(len(msg))+'s',msg)  # Message
    buf += struct.pack('B',0)  # Null termination
    try:
        #Set the whole string
        n = s.send(buf)
        print 'Number of bytes sent',n
    except:
        print 'Error sending buffer'
    try:
        s.settimeout(1.0)
        # receive data from client (data, addr)
        d = s.recv(1024)
#        reply = d[0]
#        addr = d[1]
         
        print 'Server reply : ' + d
     
    except:
        print 'Error receiving response'

    time.sleep(0.01)        
    s.close()

def get_response(msg):

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) # TCP
    s.connect((host,port))
    # Create the ethernet packet contents
    buf = rq_type['download']+rq['getresponse']
    buf += struct.pack('H',0)  # Value = 0
    buf += struct.pack('H',0)  # Index = 0
    buf += struct.pack('H',socket.htons(len(msg)+1))  # Length including null termination
    buf += struct.pack(str(len(msg))+'s',msg)  # Message
    buf += struct.pack('B',0)  # Null termination
    try:
        #Set the whole string
        n = s.send(buf)
        print 'Number of bytes sent',n
    except:
        print 'Error sending buffer'
    try:
        s.settimeout(1.0)
        # receive data from client (data, addr)
        d = s.recv(1024)
#        reply = d[0]
#        addr = d[1]
         
        print 'Server reply : ' + d[:-2]
     
    except:
        print 'Error receiving response'
        d = 'Error'

    time.sleep(0.01)        
    s.close()
    return d

def brick_flush():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) # UDP
    s.connect((host,port))

    buf = rq_type['download']+rq['flush']
    buf += struct.pack('H',0)  # Value = 0
    buf += struct.pack('H',0)  # Index = 0
    buf += struct.pack('H',0)  # Length = 0
    s.send(buf)
    s.settimeout(1.0)
    # receive data from client (data, addr)
    try:
        d = s.recv(1)
        print 'Flush returned:',d
    except:
        print 'Error (timeout?)'
    s.close()

