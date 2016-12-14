# Implements a simple tooltip string that follows the mouse
# Toggle the tooltip on/off by clicking the mouse
#
# Returns the text object, whose attributes can be changed,
# e.g. txt.set_font(...)
#
# Call the function by creating an interactive figure and axis, e.g.
#   import matplotlib.pylab as plt
#   plt.ion()
#   fig, ax = plt.subplots(1,1)
# then calling the function
#   txt = tooptip(fig,ax,callback)
# where callback is any function that takes x,y arguments
# and returns a string.  Additional arguments can be
# passed to tooltip, which will be passed on to the callback

# Example get_str function
def get_str(x,y,*args,**kwargs):
    # Return some string based on the mouse coordinates
    return 'The string is '+str(x)[:4]+' and '+str(y)[:4]+'.'

def tooltip(fig,ax,callback,*args,**kwargs):
    fig.clicked = 0
    txt = ax.text(0,0,'')
    def onmove(event):
        if event.inaxes != ax: return
        x,y = event.xdata,event.ydata
        #print 'indexes are t=%f, f=%f'%(i, j)
        txt.set_position([x,y])
        string = callback(x,y,*args,**kwargs)
        txt.set_text(string)
        fig.canvas.draw()
	
    def onclick(event):
        if fig.clicked == 0:
            fig.mid = fig.canvas.mpl_connect('motion_notify_event', onmove)
            fig.clicked = 1
        else:
            txt.set_text('')
            fig.canvas.draw()
            fig.canvas.mpl_disconnect(fig.mid)
            fig.clicked = 0

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    return txt
