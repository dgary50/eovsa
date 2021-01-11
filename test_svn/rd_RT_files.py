# cd to directory with files, change file search string
# then cut and paste into Python to generate tsys array
# for plotting
import glob
flist = glob.glob('RT20140622190*.txt')
tsys = np.zeros(len(flist)*100,dtype='float')
fn = 0
for filen in flist:
    with open(filen,'r') as f:
        lines = f.readlines() 
        for i in range(8,len(lines),10):
            # Change slice from 0 to 1, 2, or 3 to select other antennas
            tsys[fn] = float(lines[i].split()[0])
            fn += 1

