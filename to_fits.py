from casatasks import exportfits
import glob, os
import pysftp
from astropy.io import fits
# from generate_Moment_Maps import *


# sftp = pysftp.Connection('aop4.see.chalmers.se', username='leajul', password='Athena4.00')

filepath = 'DATA/R_Leo/2025-09-27/'
filenames = glob.glob(filepath+'R_Leo*fit*')
print(filenames)

for filename in filenames:
    print(filename)
    nameout = filepath +filename[len(filepath)::]+".fits"
    print(nameout)
    if not os.path.exists(nameout):
        exportfits(imagename=filename, fitsimage=nameout, dropstokes =True, velocity = True)

    # #export the center
    #
    # parameters = [datafolder, star, band, date, 0, 3, 0, [10], None]
    # ra_center, dec_center =  get_center_of_source(parameters)
    # #print(nameout)
    # #print(ra_center, dec_center)
    # #fits_image_filename = fits.util.get_testdata_filepath(filename)
    # hdul = fits.open(nameout,mode='update')  # open a FITS file
    # hdr = hdul[0].header
    # #hdr['XCENTER'] = (ra_center,'x center of a Gaussian fit of continuum image') # Add a new RA_CENTER keyword
    # #hdr['YCENTER'] = (dec_center, 'y center of a Gaussian fit of continuum image') # Add a new DEC_CENTER keyword
    #
    # hdr.append(('XCENTER', ra_center, 'x center of a Gaussian fit of continuum image'), end=True)
    # hdr.append(('YCENTER', dec_center, 'y center of a Gaussian fit of continuum image'), end=True)
    # print(hdr['XCENTER'])
    #
    #
    # #export restoring beam
    # ia = image()
    # ia.open(filename)
    # moment = ia.moments(moments=[0])
    # angle = moment.restoringbeam()["positionangle"]["value"]
    # width=moment.restoringbeam()["major"]["value"]
    # height=moment.restoringbeam()["minor"]["value"]
    #
    # hdr.append(('RBANGLE', angle, 'x center of a Gaussian fit of continuum image'), end=True)
    # hdr.append(('RBWIDTH', width, 'y center of a Gaussian fit of continuum image'), end=True)
    # hdr.append(('RBHEIGHT', height, 'x center of a Gaussian fit of continuum image'), end=True)
    #
    #
    #
    # hdul.flush() # to save
    # hdul.close()


for log_file in glob.glob("casa-*.log"):
  os.remove(log_file)
