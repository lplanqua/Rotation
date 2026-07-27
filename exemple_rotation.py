import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
import glob
import emcee, corner

plt.rcParams['font.size'] = 16


#HARDCODE VARIABLES

nbpix = 70 #lim*2 - 1#50 # nbpix
lim = 35 # size of the image in mas
pix_size =2*(lim)/nbpix
Rstar = 30

"""
Some useful functions
"""

def model(params, pix_arr, data): #x t y
    model = velocity(pix_arr,params)
    return data - model

def log_prior(params):
    axis_radian, vsini, vsys, vexp, Rshell, alpha = params
    if vsini <0:
        return -np.inf
    if vsini >50:
        return -np.inf
    if axis_radian< 0:
        return -np.inf
    if axis_radian > 2*np.pi:
        return -np.inf
    if vexp >100:       #slow wind
        return -np.inf
    if vexp <-100:
        return -np.inf
    if Rshell <Rstar:
        return -np.inf
    if Rshell >37:
        return -np.inf
    if alpha < 0:
        return -np.inf
    return 0

# Vraisemblance (log-likelihood)
def log_likelihood(params, pix_arr, data):
    axis_radian, vsini, vsys, vexp, Rshell, alpha = params
    sigma = 0.1*np.ones(np.shape(pix_arr))
    model = velocity(pix_arr,params)
    return -0.5 * np.sum(((data - model) / sigma) ** 2 + np.log(2 * np.pi * sigma ** 2))

# Posterior = prior + likelihood
def log_posterior(params, pix_arr, data):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, pix_arr, data)



def plot_surface(vel_map, params, cmap = 'jet', boolShow =True, boolSave = True, filename = None):


    axis_radian, vsini, vsys, vexp, Rshell, alpha = params
    if filename == None:
        filename = 'Figures/' +'axis' + str(int(axis_radian*180/np.pi))+ '_vsini'+"{:.1f}".format(vsini)+'_vexp'+"{:.1f}".format(vexp)+'_vsys'+"{:.1f}".format(vsys)+'_Rshell'+"{:.1f}".format(Rshell)
        filename += '.png'

    plt.figure(figsize = (12,12))
    # h = plt.contourf(X, Y, Rs)
    plt.title(r'$\Phi$ = '+ str(int(axis_radian*180/np.pi))+r'$^\circ$, $|v_{rot}\sin(i)|$ = '+ str(np.round(vsini,1))+ r' km/s'+ '\n'+ r'$v_{exp}$ = '+ str(np.round(vexp,1))+ r' km/s, $v_{sys}$ = '+ str(np.round(vsys,1))+ r' km/s, $R_{shell}$ = '+ str(np.round(Rshell,1))+' pix')
    vel_map[vel_map== 0] = np.nan
    h = plt.imshow(vel_map, extent = [-lim, lim, -lim, lim], origin = 'lower', cmap = cmap)

    plt.axis('scaled')
    ax = plt.gca()
    ax.tick_params(direction="in", which = 'both', top = True, right = True)
    if vsini != 0:
        ax.annotate("", xytext=(1.2*Rstar*np.sin(axis_radian), -1.2*Rstar*np.cos(axis_radian)), xy=(-1.2*Rstar*np.sin(axis_radian), 1.2*Rstar*np.cos(axis_radian)),
                    arrowprops=dict(arrowstyle="->",facecolor='black', lw = 2))
    ax.minorticks_on()
    cbar = plt.colorbar(label = 'RV [km/s]', pad = 0.001,  fraction=0.01,aspect=100)
    # cbar.add_lines([vsys])
    cbar.ax.axhline(y=vsys, c='k')
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


def write_to_fits(vel_map, params, filename = None):
    axis_radian, vsini, vsys, vexp, Rshell, alpha = params
    if filename == None:
        filename = 'Figures/' + 'axis' + str(int(axis_radian*180/np.pi))+ '_vsini'+"{:.1f}".format(vsini)+'_vexp'+"{:.1f}".format(vexp)+'_vsys'+"{:.1f}".format(vsys)+'_Rshell'+"{:.1f}".format(Rshell)
        filename += '.fits'
    hdr = fits.Header()
    hdr['BUNIT'] = 'km/s'
    hdr['BTYPE'] = 'velocity'

    hdr['CRVAL1'] = 0
    hdr['CRVAL2'] = 0

    hdr['CRPIX1'] = nbpix/2
    hdr['CRPIX2'] = nbpix/2

    hdr['CRDELT1'] = -pix_size #should be negative for RA axis positive toward left
    hdr['CRDELT2'] = pix_size

    # the units
    hdr['CUNIT1'] = 'mas'
    hdr['CUNIT1'] = 'mas'

    hdu = fits.PrimaryHDU(vel_map, hdr)

    hdu.writeto(filename, overwrite=True)
    return True



def velocity(pix_arr,params):
    axis_radian, vsini, vsys, vexp, Rshell, alpha = params

    Y, X = np.unravel_index(pix_arr, (nbpix,nbpix)) #X and Y swapped to speak in cartesian coordinates
    X = X*pix_size -lim
    Y = Y*pix_size -lim

    R = (X*X + Y*Y)**0.5
    #
    mask = R > Rstar
    # X[mask] = np.nan
    # Y[mask] = np.nan

    axis_perp = X*np.cos(axis_radian) + Y*np.sin(axis_radian) #x'
    axis_para =  -X*np.sin(axis_radian) + Y*np.cos(axis_radian)#y'
    vrot = axis_perp*vsini/Rshell
    vdiff = (np.abs(axis_para)*alpha/Rshell)*np.sign(vrot)
    vtot = vrot + vsys - vexp*(Rshell - R)/Rshell #+ vdiff
    vtot[mask] = 0#np.nan

    return vtot.ravel()




tot_pix = nbpix*nbpix
pix_arr = np.arange(tot_pix)

"""
The model parameters:
"""
axis = 45 #in degrees
axis_radian = axis*np.pi/180

vsini = 4 #km/s, should be positive
vsys = 0#4.9 #km/s
vexp = 0 #km/s, positive value for outflows
Rshell = 36
# pixel
alpha = -0.1
params = [axis_radian, vsini, vsys, vexp, Rshell, alpha]


vel_map_1d = velocity(pix_arr,params)

#back to 2d-array for the visualization
vel_map = np.reshape(vel_map_1d, (nbpix, nbpix))
# print(vel_map_1d)

"""
plot the result
"""
plot_surface(vel_map, params, boolShow= True, boolSave = True, cmap = 'jet',filename = 'input.png')

"""
export the result into a fits file
"""
write_to_fits(vel_map, params, filename = 'test.fits')


"""
read fits file and extract data
"""
# DIR = 'Figures/'

# params = [45*np.pi/180, 5, 0.0, 2.0, 36.0]
# axis_radian, vsini, vsys, vexp, Rshell = params

filename = 'test.fits'

# filename = glob.glob(DIR+'*.fits')[0]
print(filename)

with fits.open(filename) as hdul:
   data = hdul[0].data
   hdr = hdul[0].header

#need to extract the parameters for the image

# print(data)
# plot_surface(data, params, cmap = 'bwr',boolShow= True, boolSave = False)

data_1d = np.ravel(data)
print(np.shape(data_1d))

x0 = np.ones(np.shape(data_1d))
x0 = [0*np.pi/180, 0, 0.0, 0, Rstar, 0]


#least_squares

res_lsq = least_squares(model, x0, args=(pix_arr, data_1d))

print(res_lsq.x)
print('PA rotation axis = ', res_lsq.x[0]*180/np.pi)

res_vel_map_1d = velocity(pix_arr,res_lsq.x)
plot_surface(np.reshape(res_vel_map_1d, (nbpix, nbpix)), res_lsq.x,  cmap = 'bwr',boolShow= True, boolSave = True, filename = 'least_squares_fit.png')
err = 0.11
chi2 = np.sum(((data_1d - res_vel_map_1d)**2)/err**2)
print('chi2 = ', chi2)



#MCMC fit
ndim = len(params)
nwalkers = 50
initial_guess = x0
pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(pix_arr, data_1d))
print("Running MCMC...")
sampler.run_mcmc(pos, 15000, progress=True)
#Results
burn_in = 1000
samples = sampler.get_chain(discard=burn_in, flat=True)


params_MCMC = np.mean(samples, axis=0)
print(f"Estimation MCMC : Phi = {params_MCMC[0]:.3f}, vsini = {params_MCMC[1]:.3f}, vsys = {params_MCMC[2]:.3f}, vexp = {params_MCMC[3]:.3f}, Rshell = {params_MCMC[4]:.3f}")

MCMC_vel_map_1d = velocity(pix_arr,params_MCMC)
plot_surface(np.reshape(MCMC_vel_map_1d, (nbpix, nbpix)), params_MCMC,  cmap = 'bwr',boolShow= True, boolSave = True, filename = 'MCMC_fit.png')


samples[:,0] = samples[:,0]*180/np.pi

fig = corner.corner(
    samples,
    labels=[r"$\Phi$", r"$|v_{rot}\sin(i)|$", r"$v_{sys}$", r"$v_{exp}$", r"$R_{shell}$"],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_fmt=".2f"
)


plt.savefig('Corner.png', bbox_inches = 'tight')
plt.show()

