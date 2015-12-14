from show_capture_ import *
import numpy as np

def img2coeff(img):
    blah = np.zeros([4096,50],'float')
    for i in range(50):
        blah[:,i] = img[:,i]/median(img[:,i])
    y = np.median(blah,axis=1)
    y64 = np.zeros(64,'float')
    dx = 64
    x64 = (np.arange(64)*dx + dx/2).astype('int')
    for i in range(64):
        j = i*dx
        y64[i] = median(y[j:j+dx+1])
    figure(3)
    plot(y)
    plot(x64,y64,'ro')
    yscale('log')
    return median(y64)/y64

def calc_coeff(file):
    buf = ''
    imp, img = show_image(file,chan=0,bid=0)
    c1x = img2coeff(img)
    for i in range(64):
        buf += struct.pack('<f',c1x[i])

    imp, img = show_image(file,chan=1,bid=0)
    c1y = img2coeff(img)
    for i in range(64):
        buf += struct.pack('<f',c1y[i])

    imp, img = show_image(file,chan=2,bid=0)
    c7x = img2coeff(img)
    for i in range(64):
        buf += struct.pack('<f',c7x[i])

    imp, img = show_image(file,chan=3,bid=0)
    c7y = img2coeff(img)
    for i in range(64):
        buf += struct.pack('<f',c7y[i])

    imp, img = show_image(file,chan=2,bid=1)
    c8x = img2coeff(img)
    for i in range(64):
        buf += struct.pack('<f',c8x[i])

    imp, img = show_image(file,chan=3,bid=1)
    c8y = img2coeff(img)
    for i in range(64):
        buf += struct.pack('<f',c8y[i])

    f = open('Coeff.dat','wb')
    f.write(buf)
    f.close()

    return c1x, c1y, c7x, c7y, c8x, c8y
