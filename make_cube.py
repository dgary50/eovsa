import numpy
import pyfits

def make_cube(innames, outname):
    data = []
    for fname in innames:
        img = pyfits.getdata(fname)
        data.append(img)
    fits = pyfits.PrimaryHDU(numpy.array(data))
    fits.writeto(outname, clobber=True)

if __name__ == "__main__":
    ''' Usage:
          python make_cube.py <match string> <outfile name>
        Example:
          python make_cube.py "/lustre/solar/dgary/images/08*image.fits" "08cube.fits" 
    '''
    import glob, sys
    if len(sys.argv) == 3:
        try:
            innames = glob.glob(sys.argv[1])
            outname = sys.argv[2]
        except:
            print 'Cannot interpret',sys.argv[1],'as a valid string.'
            exit()
    else:
        print 'Incorrect number of arguments.  Need match string and outfile name'
        exit()
    make_cube(innames,outname)
    