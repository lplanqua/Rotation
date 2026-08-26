import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
import glob
import emcee, corner
import config



"""
Some useful functions for the fit
"""



def velocity(pix_arr,params):
    axis_radian, vsini, vsys, vexp = params
    nbpix = config.nbpix
    lim =  config.lim
    pix_size =  config.pix_size
    Rstar =  config.Rstar
    Rshell =  config.Rshell

    Y, X = np.unravel_index(pix_arr, (nbpix,nbpix)) #X and Y swapped to speak in cartesian coordinates
    X = X*pix_size -lim +pix_size*0.5
    Y = Y*pix_size -lim +pix_size*0.5

    R = (X*X + Y*Y)**0.5
    # print(X)
    #
    mask = R > Rstar
    # X[mask] = np.nan
    # Y[mask] = np.nan

    axis_perp = X*np.cos(axis_radian) + Y*np.sin(axis_radian) #x'
    axis_para =  -X*np.sin(axis_radian) + Y*np.cos(axis_radian)#y'
    vrot = axis_perp*np.abs(vsini)/Rstar# Rshell
    vtot = vrot + vsys - vexp*(Rshell - R)/Rshell
    vtot[mask] = 0#np.nan
    # h = plt.imshow(np.reshape(R, (nbpix, nbpix)), extent = [-lim, lim, -lim, lim], origin = 'lower', cmap = 'bwr')
    # plt.show()
    # print(stop)

    return vtot.ravel()




def model(params, pix_arr, data): #x t y
    model = velocity(pix_arr,params)
    return data - model

def log_prior(params):
    axis_radian, vsini, vsys, vexp = params
    if vsini <0:
        return -np.inf
    if axis_radian< 0:#-0.5*np.pi:
        return -np.inf
    if axis_radian > 2.5*np.pi:
        return -np.inf
    # if vexp >100:       #slow wind
    #     return -np.inf
    # if vexp <-100:
    #     return -np.inf
    return 0

# Vraisemblance (log-likelihood)
def log_likelihood(params, pix_arr, data, err):
    axis_radian, vsini, vsys, vexp = params
    # sigma = err_scalar*np.ones(np.shape(pix_arr))
    sigma = err
    model = velocity(pix_arr,params)
    return -0.5 * np.sum(((data - model) / sigma) ** 2 + np.log(2 * np.pi * sigma ** 2))

# Posterior = prior + likelihood
def log_posterior(params, pix_arr, data, err):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, pix_arr, data, err)



"""
Some useful functions for the plots
"""


def draw_self_loop(center, radius, rota=-30):

    # Add the ring
    rwidth = 0.2
    width = radius
    height = radius/2

    ring = mpatches.Arc(center, width, height, angle=rota, theta1 = 120, theta2 = 300, lw = 3)
    # Triangle edges
    offset = 0.2
    rota_rad = rota*np.pi/180 #-90 -rota*np.pi/180
    xcent  = center[0] - np.sin(rota_rad)*height/2 # + (rwidth/2)
    ycent  = center[1] - np.cos(rota_rad)*height/2 # + (rwidth/2)

    left   = [xcent + offset*np.sin(rota_rad), ycent - offset*np.cos(rota_rad)]
    right  = [xcent - offset*np.sin(rota_rad), ycent + offset*np.cos(rota_rad)]

    bottom = [offset*np.cos(rota_rad) + (left[0]+right[0])/2., ycent-offset*np.sin(rota_rad)]
    arrow  = plt.Polygon([left, right, bottom, left])
    print(left, right, bottom)
    p = PatchCollection(
        [ring, arrow],
        edgecolor = 'lime',
        facecolor = 'None',
        lw = 5
    )
    return p

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', np.nan)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

def plot_surface(vel_map, params,core_params, cmap = 'jet', chi2r = None, boolShow =True, boolSave = True, filename = None, title = None):


    axis_radian, vsini, vsys, vexp = params
    nbpix, pix_size, lim, Rstar = core_params
    Rshell = config.Rshell
    if filename == None:
        filename = 'Figures/' +'axis' + str(int(axis_radian*180/np.pi))+ '_vsini'+"{:.1f}".format(vsini)+'_vexp'+"{:.1f}".format(vexp)+'_vsys'+"{:.1f}".format(vsys)+'_Rshell'+"{:.1f}".format(Rshell)
        filename += '.png'

    plt.figure(figsize = (12,12))
    # h = plt.contourf(X, Y, Rs)
    if title == None:
        if params[1] == 0:
            plt.title('INPUT')
        else:
            plt.title(r'$\Phi$ = '+ str(int(axis_radian*180/np.pi))+r'$^\circ$, $|v_{rot}\sin(i)|$ = '+ str(np.round(vsini,1))+ r' km/s'+ '\n'+ r'$v_{exp}$ = '+ str(np.round(vexp,1))+ r' km/s, $v_{sys}$ = '+ str(np.round(vsys,1))+ r' km/s, $R_{shell}$ = '+ str(np.round(Rshell,1))+' pix')
    else:
        plt.title(title)
    # vel_map= np.array(vel_map, dtype=float)
    Cmap = plt.get_cmap(cmap)
    Cmap.set_bad(color='gray', alpha=0.25)

    vel_map[vel_map== 0] = np.nan
    #padding
    vel_map_pad = np.pad(vel_map, 1, pad_with)
    vcenter =vsys
    if config.starname == 'Betelgeuse':
        vcenter = 0
    h = plt.imshow(vel_map_pad, extent = [-lim- pix_size, +lim+pix_size, -lim-pix_size, lim+pix_size], origin = 'lower', cmap = Cmap, norm=colors.CenteredNorm(vcenter = vcenter))

    plt.axis('scaled')
    ax = plt.gca()

    ax.tick_params(direction="in", which = 'both', top = True, right = True)
    if vsini != 0:
        # print(1.2*Rstar*np.sin(axis_radian), -1.2*Rstar*np.cos(axis_radian), -1.2*Rstar*np.sin(axis_radian), 1.2*Rstar*np.cos(axis_radian))
        ax.annotate("", xytext=(1.1*Rstar*np.sin(axis_radian), -1.1*Rstar*np.cos(axis_radian)), xy=(-1.1*Rstar*np.sin(axis_radian), 1.1*Rstar*np.cos(axis_radian)),
                    arrowprops=dict(arrowstyle="->",facecolor='black', lw = 3))
        p = draw_self_loop(center=(-1.05*Rstar*np.sin(axis_radian), 1.05*Rstar*np.cos(axis_radian)),radius=2, rota = axis_radian*180/np.pi)
        # ax.add_collection(p)
        if chi2r != None:
            ax.annotate(r"$\chi^2_r = $"+"{:.2f}".format(chi2r),xy = (0.8,0.05), xycoords = 'axes fraction')
    ax.minorticks_on()
    cbar = plt.colorbar(label = 'RV [km/s]', pad = 0.001,  fraction=0.01,aspect=100)
    # cbar.add_lines([vsys])
    cbar.ax.axhline(y=vsys, c='lime', lw = 5)
    # cbar.axvline(y=vsys)
    plt.xlabel(r'$\Delta \alpha$ [mas] (<--E)')
    plt.ylabel(r'$\Delta \delta$ [mas] (N-->)')
    if boolSave:
        plt.savefig(filename, bbox_inches = 'tight')
        print('saved as : ',filename)
    if boolShow:
        plt.show()
    else:
        plt.close()
    return True


def write_to_fits(vel_map, params, core_params, filename = None):
    axis_radian, vsini, vsys, vexp = params
    nbpix, pix_size, lim, Rstar = core_params
    if filename == None:
        filename = 'Figures/' + 'axis' + str(int(axis_radian*180/np.pi))+ '_vsini'+"{:.1f}".format(vsini)+'_vexp'+"{:.1f}".format(vexp)+'_vsys'+"{:.1f}".format(vsys)+'_Rshell'+"{:.1f}".format(Rshell)
        filename += '.fits'
    hdr = fits.Header()
    hdr['BUNIT'] = 'km/s'
    hdr['BTYPE'] = 'velocity'

    hdr['CRVAL1'] = 0 #nbpix/2
    hdr['CRVAL2'] = 0 #nbpix/2

    hdr['CRPIX1'] = nbpix/2
    hdr['CRPIX2'] = nbpix/2

    hdr['CDELT1'] = pix_size #should be negative for RA axis positive toward left
    hdr['CDELT2'] = pix_size

    # the units
    hdr['CUNIT1'] = 'mas'
    hdr['CUNIT1'] = 'mas'

    hdu = fits.PrimaryHDU(vel_map, hdr)

    hdu.writeto(filename, overwrite=True)
    return True


def read_fits(filename, Rstar):
    print('Using file ',filename)
    with fits.open(filename) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header
    #if stokes params
    if data.ndim == 4:
        data = data[0][0]
    elif data.ndim == 3:
        data = data[0]
    data = np.nan_to_num(data) # Replace NaN values with 0
    nbpix = np.shape(data)[0]
    pix_size = 3
    #pix_size =2*(lim)/nbpix
    lim = nbpix*pix_size/2
    # Rstar = 21#22.6
    core_params = nbpix, pix_size, lim, Rstar

    return data, core_params


