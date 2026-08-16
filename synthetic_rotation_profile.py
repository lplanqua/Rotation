import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
import glob, os
import emcee, corner
from functions import *




#HARDCODE VARIABLES

nbpix = 30 #lim*2 - 1#50 # nbpix
lim = 31 # size of the image in mas
pix_size =2*(lim)/nbpix
Rstar = 30

"""
Some useful functions
"""

def model_dif(params, pix_arr, data): #x t y
    model = velocity_dif(pix_arr,params)
    return data - model

def log_prior_diff(params):
    axis_radian, vsini, vsys, vexp, alpha = params
    if vsini <0:
        return -np.inf
    if vsini >50:
        return -np.inf
    if axis_radian< 0:
        return -np.inf
    if axis_radian > 2*np.pi:
        return -np.inf


    return 0

# Vraisemblance (log-likelihood)
def log_likelihood_diff(params, pix_arr, data):
    axis_radian, vsini, vsys, vexp, alpha = params
    sigma = 0.1*np.ones(np.shape(pix_arr))
    model = velocity_dif(pix_arr,params)
    return -0.5 * np.sum(((data - model) / sigma) ** 2 + np.log(2 * np.pi * sigma ** 2))

# # Posterior = prior + likelihood
def log_posterior_diff(params, pix_arr, data):
    lp = log_prior_diff(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_diff(params, pix_arr, data)




def velocity_dif(pix_arr,params):
    axis_radian, vsini, vsys, vexp, alpha = params

    # Y, X = np.unravel_index(pix_arr, (nbpix,nbpix)) #X and Y swapped to speak in cartesian coordinates
    # X = X*pix_size -lim
    # Y = Y*pix_size -lim


    Y, X = np.unravel_index(pix_arr, (nbpix,nbpix)) #X and Y swapped to speak in cartesian coordinates
    X = X*pix_size -lim +pix_size*0.5
    Y = Y*pix_size -lim +pix_size*0.5

    R = (X*X + Y*Y)**0.5
    #
    mask = R > Rstar
    # X[mask] = np.nan
    # Y[mask] = np.nan

    axis_perp = X*np.cos(axis_radian) + Y*np.sin(axis_radian) #x'
    axis_para =  -X*np.sin(axis_radian) + Y*np.cos(axis_radian)#y'
    vrot = axis_perp*vsini/Rshell
    azim = np.atan(axis_para/Rstar)



    # vrot = vrot*(-1+2*np.heaviside(alpha,1) + alpha*np.sin(azim)**2)

    if alpha == 0:
        vrot = vrot
    elif alpha > 0:
        vrot = vrot*(1+ (alpha-1)*np.sin(azim)**2)/alpha
    elif alpha < 0:
        vrot = vrot*(1+(-alpha-1)*(1-np.sin(azim)**2))/(-alpha)
    # print(vrot)

    vrot += 1E-32
    vtot = vrot + vsys - vexp*(Rshell - R)/Rshell #+ vdiff
    # print(stop)
    vtot[mask] = 0#np.nan

    return vtot.ravel()




DIR = 'Figures/vel_diff_model/'

if not os.path.exists(DIR):
    os.makedirs(DIR)

tot_pix = nbpix*nbpix
pix_arr = np.arange(tot_pix)


core_params = nbpix, pix_size, lim, Rstar

"""
The model parameters:
"""
axis = 50 #in degrees
axis_radian = axis*np.pi/180

vsini = 4 #km/s, should be positive
vsys = 0#4.9 #km/´s
vexp = 0 #km/s, positive value for outflows
Rshell = 30
# pixel
alpha = -4
params = [axis_radian, vsini, vsys, vexp, alpha]


vel_map_1d = velocity_dif(pix_arr,params)

#back to 2d-array for the visualization
vel_map = np.reshape(vel_map_1d, (nbpix, nbpix))


"""
plot the result
"""
title = r'$\Phi$ = '+ str(int(axis_radian*180/np.pi))+r'$^\circ$, $|v_{rot}\sin(i)|$ = '+ str(np.round(vsini,1))+ r' km/s'+ '\n'+ r'$v_{exp}$ = '+ str(np.round(vexp,1))+ r' km/s, $v_{sys}$ = '+ str(np.round(vsys,1))+ r' km/s, $\alpha$ = '+ str(np.round(alpha,1))
filename = DIR + f"alpha = {alpha:.1f}.png" #'input.png'
plot_surface(vel_map, params[0:4],core_params, boolShow= True, boolSave = True, cmap = 'bwr',title = title, filename = filename)

"""
export the result into a fits file
"""
write_to_fits(vel_map, params[0:4],core_params, filename = DIR + 'test.fits')


"""
read fits file and extract data
"""
# DIR = 'Figures/'

# params = [45*np.pi/180, 5, 0.0, 2.0, 36.0]
# axis_radian, vsini, vsys, vexp, Rshell = params

filename = DIR +'test.fits'

# filename = glob.glob(DIR+'*.fits')[0]
print(filename)

with fits.open(filename) as hdul:
   data = hdul[0].data
   hdr = hdul[0].header

#need to extract the parameters for the image


# plot_surface(data, params, cmap = 'bwr',boolShow= True, boolSave = False)

data_1d = np.ravel(data)
data_1d = np.nan_to_num(data_1d)
# data_1d[data_1d == 'nan'] = 0



x0 = np.ones(np.shape(data_1d))
x0 = [10*np.pi/180, 1, 1, 1, 0]


#least_squares

res_lsq = least_squares(model_dif, x0, args=(pix_arr, data_1d))

print(res_lsq.x)
print('PA rotation axis = ', res_lsq.x[0]*180/np.pi)

res_vel_map_1d = velocity_dif(pix_arr,res_lsq.x)
plot_surface(np.reshape(res_vel_map_1d, (nbpix, nbpix)), res_lsq.x[0:4], core_params, cmap = 'bwr',boolShow= True, boolSave = True, filename = DIR +'least_squares_fit.png')
err = 0.11
chi2 = np.sum(((data_1d - res_vel_map_1d)**2)/err**2)
print('chi2 = ', chi2)



#MCMC fit
ndim = len(params)
nwalkers = 50
initial_guess = x0
pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_diff, args=(pix_arr, data_1d))
print("Running MCMC...")
sampler.run_mcmc(pos, 15000, progress=True)
#Results
burn_in = 1000
samples = sampler.get_chain(discard=burn_in, flat=True)


params_MCMC = np.mean(samples, axis=0)
print(f"Estimation MCMC : Phi = {params_MCMC[0]:.3f}, vsini = {params_MCMC[1]:.3f}, vsys = {params_MCMC[2]:.3f}, vexp = {params_MCMC[3]:.3f}, alpha = {params_MCMC[4]:.3f}")

MCMC_vel_map_1d = velocity_dif(pix_arr,params_MCMC)
plot_surface(np.reshape(MCMC_vel_map_1d, (nbpix, nbpix)), params_MCMC[0:4], core_params, cmap = 'bwr',boolShow= True, boolSave = True, filename = DIR +'MCMC_fit.png')


samples[:,0] = samples[:,0]*180/np.pi

fig = corner.corner(
    samples,
    labels=[r"$\Phi$", r"$|v_{rot}\sin(i)|$", r"$v_{sys}$", r"$v_{exp}$", r"$\alpha$"],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_fmt=".2f"
)


plt.savefig(DIR +'Corner.png', bbox_inches = 'tight')
plt.show()

