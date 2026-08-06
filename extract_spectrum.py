import os, glob
from pathlib import Path
from time import time
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Ellipse
from astropy.wcs import WCS
import matplotlib.colors as colors
import astropy.units as u
from astropy.io import fits
from astropy.modeling import models, fitting
import config




from matplotlib.gridspec import GridSpec

plt.rcParams['font.size'] = 16



starname = 'R_Leo' #or R_Dor, R_Leo or Betelgeuse
date = '2025-07-15'



if starname == 'Betelgeuse':
    nbpix = 15 #lim*2 - 1#50 # nbpix
    lim = 22.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    Rstar = 21#22.6
    Rshell = 40
    vsys0 = 0

elif starname == 'R_Dor':
    # nbpix = 17 #lim*2 - 1#50 # nbpix
    # lim = 25.5 # size of the image in mas
    # pix_size =3 #*(lim)/nbpix
    Rstar = 24#22.6
    Rshell = 30
    vsys0 = 9.5


elif starname == 'R_Leo':
    nbpix = 15 #lim*2 - 1#50 # nbpix
    lim = 22.5 # size of the image in mas
    pix_size =2*(lim)/nbpix

    Rstar = 21#22.6
    Rshell = 30
    vsys0 = 10



if starname == 'Betelgeuse':
    transition = '28SiOv2'
    # transition = '12CO'
    date = ''
    input_dir =  'DATA/'+starname+ '/'

elif starname == 'R_Dor' or starname == 'R_Leo':
    transition = 'SiO_v=2_8-7'
    transition = 'SO_3_v=1_8-7'
    input_dir =  'DATA/'+starname+ '/' + date + '/'




def extract_aperture(maps, px,py, aperture):
    XMESH, YMESH = np.meshgrid(np.linspace(nbpix%2, np.shape(maps)[1]-nbpix%2, np.shape(maps)[1]),np.linspace(nbpix%2, np.shape(maps)[1]-nbpix%2, np.shape(maps)[1]))
    mask = ((XMESH - px)**2 + (YMESH - py)**2 < (aperture)**2)
    mask = mask[ np.newaxis,:, :]

    data = np.where(mask, maps, np.nan)
    # spectrum = np.nansum(data, axis = (1,2))*pix*pix/(rb_height*rb_width*1.133)#maps[:,py,px]
    # if starname == 'Betelgeuse':
    spectrum = np.nansum(data, axis = (1,2))#maps[:,py,px]
    data = np.sum(data[:,:,:], axis = 0)

    return spectrum, data






#exracting the cube from the model

filename = glob.glob(input_dir + '*'+ transition+'*model.image*.fits')[0]

residual = glob.glob(input_dir + '*'+ transition+'*model.residual*.fits')[0]
res = fits.open(residual)[0].data
print(filename)

rb_height= rb_width = 1
ia = fits.open(filename)
hdr = ia[0].header


# rb = ia[0].data
if starname == 'Betelgeuse':
    rb_height = 1
    rb_width = 1
else:
    rb_height = hdr['BMAJ']*1E3
    rb_width = hdr['BMIN']*1E3


maps = ia[0].data


vmax = np.max(maps)
vmin = np.min(maps)

nb_channels = np.shape(maps)[0]
print(nb_channels)

rv_min = hdr['CRVAL3']*1E-3
width = hdr['CDELT3']*1E-3
rv_max = rv_min + nb_channels*width

# print(hdr['HISTORY'])


channels = np.arange(rv_min, rv_max,width)


# arg_channels = np.argwhere((channels >=-20) & (channels <= 20))



nbpix =  np.shape(maps)[1]


x0 = np.shape(maps)[1]/2
y0 = np.shape(maps)[1]/2
center = x0
print(x0, y0)

leftlim = x0
rightlim = x0-nbpix
bottomlim = -y0
toplim = nbpix-y0

pix = config.pix_size
leftlim *= pix
rightlim *=pix
bottomlim *=pix
toplim *=pix


#info on the aperture
delta_xB = +5
delta_yB = 0

delta_xA = -6
delta_yA = 3

apertureB = 0.5 # in px
apertureA = 0.5
# aperture = # in px
px = center-delta_xB #columns
py = center-delta_yB #rows
# print(px, py)

spectrumB, dataB = extract_aperture(maps, center-delta_xB,center-delta_yB, apertureB )
spectrumA, dataA = extract_aperture(maps,  center-delta_xA,center-delta_yA, apertureA )



res_spectrumB, _ = extract_aperture(res, center-delta_xB,center-delta_yB, apertureB )
res_spectrumA, _ = extract_aperture(res,  center-delta_xA,center-delta_yA, apertureA )


# fig = plt.figure(figsize = (12,40))
fig, axs = plt.subplots(1,2,figsize=(14, 6), constrained_layout = True)
# ax = plt.gca()
cmap = 'hot'
Cmap = plt.get_cmap(cmap)
Cmap.set_bad(color='gray', alpha=0.25)
img = axs[0].imshow(np.sum(maps[:,:,:], axis = 0), cmap = Cmap,origin = 'lower', extent = [leftlim, rightlim, bottomlim, toplim])
fig.colorbar(img, ax=axs[0], pad = 0.001,  fraction=0.01,aspect=100)
#data = data[20,:,:]
center = 0
axs[0].imshow(dataB, cmap = 'Blues_r',  origin = 'lower', alpha = 1, extent = [leftlim, rightlim, bottomlim, toplim])#, extent=extent, alpha = 0.5)
blobB = Circle(((center+delta_xB)*pix, (center-delta_yB)*pix), radius = apertureB*pix, color='blue', fill=False,linewidth=3)
axs[0].add_patch(blobB)
axs[0].imshow(dataA, cmap = 'Greens_r',  origin = 'lower', alpha = 1, extent = [leftlim, rightlim, bottomlim, toplim])#, extent=extent, alpha = 0.5)
blobA = Circle(((center+delta_xA)*pix, (center-delta_yA)*pix), radius = apertureA*pix, color='green', fill=False,linewidth=3)
axs[0].add_patch(blobA)
#axs[0].plot(py,px, '*', lw = 8, color = 'lime')


# ax = plt.gca()
axs[0].tick_params(direction="in", which = 'both', top = True, right = True, color = 'black')
axs[0].minorticks_on()
axs[0].set_xlabel(r' $\Delta$RA [mas]')
axs[0].set_ylabel(r'$\Delta$DEC [mas]')


# plt.subplot(1,2,2)


axs[1].plot(channels,spectrumB+ res_spectrumB,  color = 'blue',drawstyle='steps-mid', label = 'Blob B')
axs[1].plot(channels,spectrumB,  color = 'darkblue', ls = '--', lw = 2,label = 'Fit')
# axs[1].axhline(y = 0, ls='--', color = 'blue', alpha = 0.5)
plt.legend(loc = 2)
# ax2 = axs[1].twinx()
# secax.set_ylabel('Secondary-Y-Axis')

axs[1].plot(channels,res_spectrumA+ spectrumA,  color = 'green',drawstyle='steps-mid', label = 'Blob G')
axs[1].plot(channels,spectrumA,  color = 'darkgreen', ls = '--', lw = 2,label = 'Fit')
#plt.plot(xfit, fit_y, '-', color = 'red')
plt.legend(loc = 3)
axs[1].axvline(x = 0, ls='--', color = 'black', alpha = 0.5)
axs[1].axhline(y = 0, ls='--', color = 'black', alpha = 0.5)
axs[1].set_ylabel(r'Flux density [Jy/beam]')
axs[1].set_xlabel(r'RV [km/s]')
axs[1].tick_params(direction="in", which = 'both', top = True, right = True)
axs[1].minorticks_on()
# axs[1].yaxis.set_ticklabels([])



plt.savefig(input_dir + starname+ '_'+ transition + '_' + date + '_extract_spectrum.png')
plt.show()




# for log_file in glob.glob("casa-*.log"):
#     os.remove(log_file)

