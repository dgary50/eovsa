import get_sat_info as gs
import pcapture2 as p

def xy2rl(filename, satname):
    out = p.rd_jspec(filename)
    sat = gs.get_sat_info([satname])[0]
    freq = np.linspace(0,4095,4096)*0.4/4096 + 12.15
    frq = sat['freqlist']
    pol = sat['pollist']
    freqmhz = (freq*10000. + 0.5).astype('int')/10.
    ridx, = np.where(pol == 'R')
    lidx, = np.where(pol == 'L')
    rfrq = frq[ridx]
    ridx = []
    for f in rfrq:
        try:
            ridx.append(np.where(f == freqmhz)[0][0])
        except:
            pass
    ridx = np.array(ridx)

    nant = 13
    ipol = np.zeros((nant,4096),dtype='float')
    vpol = np.zeros((nant,4096),dtype='float')
    rpol = np.zeros((nant,4096),dtype='float')
    lpol = np.zeros((nant,4096),dtype='float')
    for k in range(nant):
        xx,yy,xy,yx = out['a'][k,:,:,30]

        pfitr = np.polyfit(freq[ridx],unwrap(angle(xy[ridx])),1)
        phi = polyval(pfitr,freq)
        xyp = xy*(cos(phi+pi/2)-1j*sin(phi+pi/2))
        yxp = yx*(cos(phi+pi/2)+1j*sin(phi+pi/2))

        ipol[k] = abs(xx+yy)
        vpol[k] = imag(yxp - xyp)
        rpol[k] = real(xx+yy) + imag(yxp - xyp)
        lpol[k] = real(xx+yy) - imag(yxp - xyp)

    # Normalize to antenna 6 IPOL (kind of random)
    fac = ipol/ipol[5,:]
    f, ax = subplots(4,1)
    for k in range(nant):
        ax[0].plot(freq,ipol[k]/fac[k])
        ax[1].plot(freq,vpol[k]/fac[k])
        ax[2].plot(freq,rpol[k]/fac[k])
        ax[3].plot(freq,lpol[k]/fac[k])

    ax[0].text(0.05,0.8,'Stokes I',transform=ax[0].transAxes,fontsize=12)
    ax[1].text(0.05,0.8,'Stokes V',transform=ax[1].transAxes,fontsize=12)
    ax[2].text(0.05,0.8,'RCP',transform=ax[2].transAxes,fontsize=12)
    ax[3].text(0.05,0.8,'LCP',transform=ax[3].transAxes,fontsize=12)

    for i in range(4):
        yrng = ax[i].yaxis.get_data_interval()
        ax[i].set_xlim(ax[i].xaxis.get_data_interval())
        for j in range(len(frq)):
            if pol[j] == 'R':
                ax[i].plot(frq[j]*ones(2)/1000.,yrng,color='red',linewidth=2)
            else:
                ax[i].plot(frq[j]*ones(2)/1000.,yrng,color='green',linewidth=2)
    f.suptitle(sat['name']+' Communication Satellite',fontsize=18)
    ax[3].set_xlabel('Frequency [GHz]')
    if sat['name'] == 'Ciel-2':
        ax[3].annotate('Beacon',(12.209,14000),xytext=(12.17,30000),arrowprops=dict(width=2,headwidth=6,frac=0.2,facecolor='black', shrink=0.05))
    show()
    return ipol,vpol,rpol,lpol
