from casatasks import exportfits
import glob, os
# from astropy.io import fits
# from generate_Moment_Maps import *

starname = 'R_Dor'
date = '2023-10-12'
filepath = 'DATA/'+ starname+ '/' + date + '/'
filenames = glob.glob(filepath+ starname+'*fit*')
print(filenames)

for filename in filenames:
    # print(filename)
    if '.fits' in filename:
      continue
    nameout = filepath +filename[len(filepath)::]+".fits"
    print(nameout)
    if not os.path.exists(nameout):
        exportfits(imagename=filename, fitsimage=nameout, dropstokes =True, velocity = True)



for log_file in glob.glob("casa-*.log"):
  os.remove(log_file)
