import numpy as np
from functions import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
from astropy.io import fits
import glob

plt.rcParams['font.size'] = 16

starname = 'R_Leo' # R_Dor , R_Leo or Betelgeuse
date = '2025-07-15'

#HARDCODE VARIABLES

if starname == 'Betelgeuse':
    nbpix = 15 #lim*2 - 1#50 # nbpix
    lim = 22.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    print(pix_size)
    Rstar = 21#22.6
    Rshell = 40
    vsys0 = 0

elif starname == 'R_Dor':
    nbpix = 17 #lim*2 - 1#50 # nbpix
    lim = 25.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    print(pix_size)
    Rstar = 24#22.6
    Rshell = 30
    vsys0 = 9.5


elif starname == 'R_Leo':
    nbpix = 15 #lim*2 - 1#50 # nbpix
    lim = 22.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    print(pix_size)
    Rstar = 21#22.6
    Rshell = 30
    vsys0 = 10



if starname == 'Betelgeuse':
    transition = '28SiOv2'
    # transition = '12CO'
    input_dir =  'DATA/'+starname+ '/'

elif starname == 'R_Dor' or starname == 'R_Leo':
    transition = 'SiO_v=3_8-7.clean.large_scale' #'SO_3_v=1_8-7' 'SiO_v=2_8-7'
    transition = 'SiO_v=2_8-7'
    input_dir =  'DATA/'+starname+ '/' + date + '/'




amp_file = glob.glob(input_dir + '*'+ transition+'*.amp_*.fits')[0]
vel_file = glob.glob(input_dir + '*'+ transition+'*center.vel*.fits')[0]
fwhm_file = glob.glob(input_dir + '*'+ transition+'*.fwhm_*.fits')[0]

amp_file_err = glob.glob(input_dir + '*'+ transition+'*.amp.err*.fits')[0]
vel_file_err = glob.glob(input_dir + '*'+ transition+'*center.err*.fits')[0]
fwhm_file_err = glob.glob(input_dir + '*'+ transition+'*.fwhm.err*.fits')[0]



print(amp_file, amp_file_err, vel_file, vel_file_err, fwhm_file, fwhm_file_err)
amp_data, core_params = read_fits(amp_file, Rstar)
vel_data, _ = read_fits(vel_file, Rstar)
fwhm_data, _ = read_fits(fwhm_file, Rstar)

amp_err_data, _ = read_fits(amp_file_err, Rstar)
vel_err_data, _ = read_fits(vel_file_err, Rstar)
fwhm_err_data, _ = read_fits(fwhm_file_err, Rstar)


def create_subplot(ax, data, title, cmap):
    ax.set_title(title)
    Cmap = plt.get_cmap(cmap)
    Cmap.set_bad(color='gray', alpha=0.25)

    data[data== 0] = np.nan
    img = ax.imshow(data, origin = 'lower', cmap = Cmap)
    # ax1.colorbar(label = 'RV [km/s]', pad = 0.001,  fraction=0.01,aspect=100)
    ax.tick_params(direction="in", which = 'both', top = True, right = True)
    ax.minorticks_on()
    fig.colorbar(img, ax=ax, pad = 0.001,  fraction=0.01,aspect=100)
    stellarR = Circle((np.shape(data)[0]/2, np.shape(data)[0]/2), radius = 7.5, color='lime', fill=False,linewidth=2)
    ax.add_patch(stellarR)
    ax.set_xticklabels([])
    ax.set_yticklabels([])


    return True


fig = plt.figure(figsize = (14,8), layout="constrained")

gs = gridspec.GridSpec(2, 3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1],sharey=ax1,  sharex=ax1)
ax3 = fig.add_subplot(gs[0, 2],sharey=ax1,  sharex=ax1)

ax4 = fig.add_subplot(gs[1, 0],sharey=ax1,  sharex=ax1)
ax5 = fig.add_subplot(gs[1, 1],sharey=ax1,  sharex=ax1)
ax6 = fig.add_subplot(gs[1, 2],sharey=ax1,  sharex=ax1)

create_subplot(ax1, amp_data, 'Amplitude [Jy/beam]', 'hot')
create_subplot(ax2, vel_data, 'Velocity [km/s]', 'bwr')
create_subplot(ax3, fwhm_data, 'FWHM [km/s]', 'jet')
create_subplot(ax4, amp_err_data, 'err(Amplitude)', 'hot')
create_subplot(ax5, vel_err_data, 'err(Velocity)', 'bwr')
create_subplot(ax6, fwhm_err_data, 'err(FWHM)', 'jet')

plt.savefig(input_dir + transition+'specfit_output.png', bbox_inches = 'tight')
plt.show()

# gs = Gridspe











